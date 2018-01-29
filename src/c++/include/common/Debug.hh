/**
 ** Isaac Genome Alignment Software
 ** Copyright (c) 2010-2017 Illumina, Inc.
 ** All rights reserved.
 **
 ** This software is provided under the terms and conditions of the
 ** GNU GENERAL PUBLIC LICENSE Version 3
 **
 ** You should have received a copy of the GNU GENERAL PUBLIC LICENSE Version 3
 ** along with this program. If not, see
 ** <https://github.com/illumina/licenses/>.
 **
 ** \file Debug.hh
 **
 ** \brief Various debugging-related helpers
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_LOG_THREAD_TIMESTAMP_HH
#define iSAAC_LOG_THREAD_TIMESTAMP_HH

#include "config.h"

#include <atomic>
#include <memory>
#include <typeinfo>

#include <boost/algorithm/string.hpp>
#include <boost/date_time.hpp>
#include <boost/io/ios_state.hpp>
#include <boost/thread.hpp>

#include "common/SystemCompatibility.hh"

namespace isaac
{
namespace common
{

#define iSAAC_PROFILING_NOINLINE
//#define iSAAC_PROFILING_NOINLINE __attribute__((noinline))


static std::ostream nostream(0);
//TODO: check why simple CerrLocker(std::cerr) << ... << is not good enough
/**
 * \brief helper macro to simplify the thread-guarded logging. All elements on a single << line are serialized
 * under one CerrLocker
 */
#define ISAAC_THREAD_CERR \
    if(const ::isaac::common::detail::CerrLocker &isaac_cerr_lock = ::isaac::common::detail::CerrLocker()); \
    else (isaac_cerr_lock.cerrBlocked() ? ::isaac::common::nostream : std::cerr) << isaac::common::detail::ThreadTimestamp()

#define ISAAC_ASSERT_CERR \
    if(::isaac::common::detail::CerrLocker isaac_cerr_lock = ::isaac::common::detail::CerrLocker()); \
    else std::cerr << isaac::common::detail::ThreadTimestamp()


#define ISAAC_SCOPE_BLOCK_CERR if (const ::isaac::common::detail::CerrBlocker blocker = ::isaac::common::detail::CerrBlocker()); else

/**
 * \brief Evaluates expression always (even if NDEBUG is set and so on). Also uses ostream serialization which,
 *        unlike the standard assert, has shown not to allocate the dynamic memory at the time when you least
 *        expect this to happen.
 */


#define ISAAC_VERIFY_MSG(expr, msg) {if (expr) {} else \
{ ISAAC_ASSERT_CERR << "ERROR: ***** Internal Program Error - assertion (" << #expr << ") failed in " \
    << (BOOST_CURRENT_FUNCTION) << ":" << __FILE__ << '(' << __LINE__ << "): " << msg << std::endl; \
    ::isaac::common::terminateWithCoreDump();}}

#if (ISAAC_BUILD_TYPE == ISAAC_BUILD_Release)
#define ISAAC_ASSERT_MSG(expr, msg)
#else
#define ISAAC_ASSERT_MSG(expr, msg) {if (expr) {} else \
{ ISAAC_ASSERT_CERR << "ERROR: ***** Internal Program Error - assertion (" << #expr << ") failed in " \
    << (BOOST_CURRENT_FUNCTION) << ":" << __FILE__ << '(' << __LINE__ << "): " << msg << std::endl; \
    ::isaac::common::terminateWithCoreDump();}}
#endif

inline std::string parseStat(const std::string &stat)
{
    std::vector<std::string> statFields;
    boost::algorithm::split(statFields, stat,  boost::algorithm::is_any_of(" "));
    return std::string(statFields.at(22) + "vm " + statFields.at(23) + "res");
}

class ScopedMallocBlock : boost::noncopyable
{
public:
    enum Mode
    {
        Invalid = 0,
        Off,
        Warning,
        Strict
    };

    ScopedMallocBlock(const Mode mode);
    ~ScopedMallocBlock();
private:
    const Mode mode_;

    friend class ScopedMallocBlockUnblock;
    void block();
    void unblock();
};

class ScopedMallocBlockUnblock : boost::noncopyable
{
    ScopedMallocBlock &block_;
public:
    ScopedMallocBlockUnblock(ScopedMallocBlock & block);
    ~ScopedMallocBlockUnblock();
};

namespace detail
{
class ThreadTimestamp
{
public:
};

struct IndentBase
{
protected:
    static iSAAC_THREAD_LOCAL unsigned width;
public:
    static unsigned getWidth() {return width;}
};

template <int i>
struct IndentT: public IndentBase
{
    IndentT() {width += i;}
    ~IndentT() {width -= i;}
};

class IndentStorer
{
public:
     IndentStorer(const IndentBase& indent) : indent(indent) {}
     const IndentBase & indent;
};

inline std::ostream & operator << (std::ostream &os, const IndentStorer & indenter) {
    os << std::setfill(' ') << std::setw(indenter.indent.getWidth()) << "";
    return os;
}

/**
 * \brief formats time stamp and thread id to simplify threaded logging
 */
inline std::ostream & operator << (std::ostream &os, const ThreadTimestamp &) {

    // IMPORTANT: this is the way to serialize date without causing any dynamic memory operations to occur
    ::std::time_t t;
    ::std::time(&t);
    ::std::tm curr, *curr_ptr;
    curr_ptr = boost::date_time::c_time::localtime(&t, &curr);

    os << (curr_ptr->tm_year + 1900) << '-' <<
        std::setfill('0') << std::setw(2) << (curr_ptr->tm_mon + 1) << '-'  <<
        std::setfill('0') << std::setw(2) << curr_ptr->tm_mday << ' '  <<

        std::setfill('0') << std::setw(2) << curr_ptr->tm_hour << ':' <<
        std::setfill('0') << std::setw(2) << curr_ptr->tm_min << ':' <<
        std::setfill('0') << std::setw(2) << curr_ptr->tm_sec << ' ' <<
        "\t[" << boost::this_thread::get_id() << "]\t";
    return os;
}

/**
 * \brief Blocks ISAAC_THREAD_CERR messages. Use for unit tests
 */
class CerrBlocker
{
    static std::atomic_int cerrBlocked_;
public:
    CerrBlocker();

    ~CerrBlocker();

    operator bool () const {
        return false;
    }

    static bool blocked() {return cerrBlocked_;}
};

struct CerrStreamBlocker : public std::ostream
{
    const bool block_;

	CerrStreamBlocker(const bool block) : std::ostream(0), block_(block) {}

    template <typename AnyT>
    const CerrStreamBlocker& operator << (const AnyT value) const
    {
        if (!block_)
        {
            std::cerr << value;
        }

        return *this;
    }

//    template <typename AnyT>
//    friend const CerrStreamBlocker& operator << (const CerrStreamBlocker& sb, const AnyT &value)
//    {
//        if (!sb.block_)
//        {
//            std::cerr << value;
//        }
//
//        return sb;
//    }


};
/**
 * \brief Guards std::cerr for the duration of CerrLocker existance
 *        Restores any changes made to ios::base
 */
class CerrLocker
{
    // some people allocate memory from under their trace code. For example by using boost::format.
    // if memory control is on, we don't want them to be dead-locked on their own thread cerrMutex_.
    static boost::recursive_mutex cerrMutex_;
    boost::lock_guard<boost::recursive_mutex> lock_;
    boost::io::ios_base_all_saver ias_;

public:

    CerrLocker(const CerrLocker &that) : lock_(cerrMutex_), ias_(std::cerr){
    }
    CerrLocker() : lock_(cerrMutex_), ias_(std::cerr) {
    }
    operator bool () const {
        return false;
    }

    bool cerrBlocked() const {return CerrBlocker::blocked();}
};

inline CerrBlocker::CerrBlocker()
{
//    ISAAC_ASSERT_CERR << "cerr blocked" << std::endl;
    ++cerrBlocked_;
}

inline CerrBlocker::~CerrBlocker()
{
    ISAAC_ASSERT_MSG(cerrBlocked_, "Attempt to unblock more times than blocked. something is really wrong");
    --cerrBlocked_;
//    ISAAC_ASSERT_CERR << "cerr unblocked" << std::endl;
}


inline void assertion_failed_msg(char const * expr, char const * msg, char const * function,
                                 char const * file, int64_t line)
{
    ISAAC_ASSERT_CERR
    << "ERROR: ***** Internal Program Error - assertion (" << expr << ") failed in "
    << function << ":" << file << '(' << line << "): " << msg << std::endl;

    common::terminateWithCoreDump();
}

} // namespace detail

inline std::time_t time()
{
    std::time_t ret;
    ISAAC_VERIFY_MSG(-1 != ::std::time(&ret), "std::time failed, errno: " << errno << strerror(errno));
    return ret;
}

} // namespace common

//isaac namespace IsaacDebugTraceIndent gets shadowed by non-incrementing type inside the first trace. This way the indent is
//consistent within the function
typedef isaac::common::detail::IndentT<1> IsaacDebugTraceIndent;


/**
 ** \brief Provide a mechanism for detailed level of debugging
 **/
#ifdef ISAAC_THREAD_CERR_DEV_TRACE_ENABLED
    #define ISAAC_THREAD_CERR_DEV_TRACE(trace) {ISAAC_THREAD_CERR << trace << std::endl;}
    #define ISAAC_DEV_TRACE_BLOCK(block) block
    #define ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(clusterId, trace) ISAAC_THREAD_CERR_DEV_TRACE(trace)
#else
    #define ISAAC_THREAD_CERR_DEV_TRACE(blah)
    #ifdef ISAAC_THREAD_CERR_DEV_TRACE_ENABLED_CLUSTER_ID
        #define ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID_INDENTING(clusterId, trace, line) const IsaacDebugTraceIndent indent ## line; typedef isaac::common::detail::IndentBase IsaacDebugTraceIndent; isaac::common::detail::IndentStorer indentStorer ## line(indent ## line); if(ISAAC_THREAD_CERR_DEV_TRACE_ENABLED_CLUSTER_ID == (clusterId)) { ISAAC_THREAD_CERR << indentStorer ## line << trace << std::endl; }
        #define ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID_PROXY(clusterId, trace, line) ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID_INDENTING(clusterId, trace, line)
        #define ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(clusterId, trace) ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID_PROXY(clusterId, trace, __LINE__)
        #define ISAAC_DEV_TRACE_BLOCK(block) block
    #else
        #define ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(clusterId, blah)
        #define ISAAC_DEV_TRACE_BLOCK(block)
    #endif
#endif

} // namespace isaac



#endif // #ifndef iSAAC_LOG_THREAD_TIMESTAMP_HH
