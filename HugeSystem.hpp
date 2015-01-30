#ifndef _SNARKLIB_HUGE_SYSTEM_HPP_
#define _SNARKLIB_HUGE_SYSTEM_HPP_

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <functional>
#include <istream>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>
#include "Rank1DSL.hpp"

namespace snarklib {

////////////////////////////////////////////////////////////////////////////////
// append/iterate over constraint system on disk or pass-through to R1System<T>
//

template <typename T>
class HugeSystem
{
    // initialize for iterating over a constraint system on disk
    template <typename T2> friend
    std::istream& operator>> (std::istream& is, HugeSystem<T2>& a);

public:
    // appending to a constraint system and write out to disk
    HugeSystem(const std::string& filePrefix,
               const std::size_t maxSize) // set to 0 for no limit
        : m_filePrefix(filePrefix),
          m_fileCount(0),
          m_constraintsPerFile(maxSize),
          m_totalConstraints(0),
          m_minIndex(-1), // maximum possible number
          m_maxIndex(0),  // minimum possible number
          m_numCircuitInputs(0),
          m_error(false)
    {}

    // iterating over a constraint system on disk
    HugeSystem(const std::string& filePrefix)
        : HugeSystem{filePrefix, 0}
    {}

    // pass-through mode to R1System<T>
    HugeSystem()
        : HugeSystem{std::string()}
    {}

    // appending to a constraint system and write out to disk
    void clearAppend(const std::string& filePrefix,
                     const std::size_t maxSize) // set to 0 for no limit
    {
        m_filePrefix = filePrefix;
        m_fileCount = 0;
        m_constraintsPerFile = maxSize;
        m_totalConstraints = 0;
        m_minIndex = -1; // maximum possible number
        m_maxIndex = 0;  // minimum possible number
        m_numCircuitInputs = 0;
        m_r1system.clear();
        m_error = false;
    }

    // iterating over a constraint system on disk
    void clearIterate(const std::string& filePrefix) {
        clearAppend(filePrefix, 0);
    }

    // pass-through mode to R1System<T>
    void clear() {
        clearIterate(std::string());
    }

    // true if there was an error while reading or writing files
    bool operator! () const { return m_error; }

    std::size_t minIndex() const { return m_minIndex; }
    std::size_t maxIndex() const { return m_maxIndex; }
    std::size_t numCircuitInputs() const { return m_numCircuitInputs; }

    const std::string& filePrefix() const { return m_filePrefix; }
    std::size_t fileCount() const { return m_fileCount; }
    std::size_t constraintsPerFile() const { return m_constraintsPerFile; }
    std::size_t totalConstraints() const { return m_totalConstraints; }

    void addConstraint(const R1Constraint<T>& d) {
        m_r1system.addConstraint(d);

        ++m_totalConstraints;

        m_minIndex = std::min(m_minIndex, m_r1system.minIndex());
        m_maxIndex = std::max(m_maxIndex, m_r1system.maxIndex());

        if (m_r1system.size() == m_constraintsPerFile) {
            flushToFile();
        }
    }

    // finished appending, write to disk
    void finalize(const std::size_t numCircuitInputs) {
        m_numCircuitInputs = numCircuitInputs;

        if (passThroughMode()) return;

        // write last block of constraints, if any
        if (0 != m_r1system.size()) {
            flushToFile();
        }

        // index file name is the prefix without a number
        std::ofstream ofs(m_filePrefix);
        if (!ofs) {
            m_error = true; // failure
        } else {
            writeIndexFile(ofs);
        }
    }

    // begin iterating, read from disk
    bool loadIndex() {
        if (passThroughMode()) return true; // ok

        // index file name is the prefix without a number
        std::ifstream ifs(m_filePrefix);
        return (!ifs)
            ? !(m_error = true) // failure
            : readIndexFile(ifs);
    }

    void swap_AB() {
        mapLambda(
            [] (R1System<T>& S) -> bool {
                S.swap_AB();
                return true; // write back to disk
            });
    }

    bool swap_AB_if_beneficial() {
        std::vector<int>
            touchA(maxIndex() + 1, false),
            touchB(maxIndex() + 1, false);

        mapLambda(
            [&touchA, &touchB] (const R1System<T>& S) -> bool {
                for (const auto& c : S.constraints()) {
                    for (const auto& t : c.a().terms()) touchA[t.index()] = true;
                    for (const auto& t : c.b().terms()) touchB[t.index()] = true;
                }

                return false; // do not write back to disk
            });

        const auto
            nzA = std::count_if(touchA.begin(), touchA.end(), [] (int i) -> bool { return i; }),
            nzB = std::count_if(touchB.begin(), touchB.end(), [] (int i) -> bool { return i; });

        return (nzA < nzB)
            ? swap_AB(), true
            : false;
    }

    bool mapLambda(std::function<bool (const R1System<T>&)> func) const {
        return passThroughMode()
            ? func(m_r1system), true // ok
            : mapLambdaFiles(func);
    }

    bool mapLambda(std::function<bool (R1System<T>&)> func) {
        return passThroughMode()
            ? func(m_r1system), true // ok
            : mapLambdaFiles(func);
    }

private:
    bool passThroughMode() const {
        return m_filePrefix.empty();
    }

    template <typename U>
    bool mapLambdaFiles(std::function<U>& func) const {
        // read files on disk
        for (std::size_t i = 0; i < m_fileCount; ++i) {
            // consecutively numbered filenames
            std::stringstream ss;
            ss << m_filePrefix << i;
            const auto& filename = ss.str();

            // use local object so this function can be const
            R1System<T> constraintSystem;

            // read in block of constraints
            {
                std::ifstream ifs(filename);
                if (!ifs || !constraintSystem.marshal_in(ifs)) {
                    return false; // failure
                }
            }

            // apply lambda and write back out to disk if required
            if (func(constraintSystem)) {
                std::ofstream ofs(filename);
                if (!ofs) return false; // failure

                constraintSystem.marshal_out(ofs);
            }
        }

        return true; // ok
    }

    void flushToFile() {
        // consecutively numbered filenames
        std::stringstream ss;
        ss << m_filePrefix << m_fileCount;
        const auto& filename = ss.str();

        // write to file stream
        std::ofstream ofs(filename);
        if (!ofs) {
            // failure
            m_error = true;

        } else {
            m_r1system.marshal_out(ofs);

            // ok, start new block
            m_r1system.clear();
            ++m_fileCount;
        }
    }

    void writeIndexFile(std::ostream& os) const {
        os << filePrefix() << std::endl
           << fileCount() << std::endl
           << m_constraintsPerFile << std::endl
           << m_totalConstraints << std::endl
           << minIndex() << std::endl
           << maxIndex() << std::endl
           << numCircuitInputs() << std::endl;
    }

    bool readIndexFile(std::istream& is) {
        return (!(is >> m_filePrefix) ||
                !(is >> m_fileCount) ||
                !(is >> m_constraintsPerFile) ||
                !(is >> m_totalConstraints) ||
                !(is >> m_minIndex) ||
                !(is >> m_maxIndex) ||
                !(is >> m_numCircuitInputs))
            ? !(m_error = true) // failure
            : true;
    }

    // file names are concatenated prefix with count
    std::string m_filePrefix;
    std::size_t m_fileCount;

    std::size_t m_constraintsPerFile, m_totalConstraints;
    std::size_t m_minIndex, m_maxIndex, m_numCircuitInputs;

    // only used for pass-through or appending, not iterating over disk
    R1System<T> m_r1system;

    // detect any errors in I/O
    bool m_error;
};

// initialize for iterating over a constraint system on disk
template <typename T>
std::istream& operator>> (std::istream& is, HugeSystem<T>& a) {
    a.clear();
    a.readIndexFile(is);
    return is;
}

} // namespace snarklib

#endif
