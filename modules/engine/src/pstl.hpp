// Referencing https://github.com/ParRes/Kernels/blob/master/Cxx11/prk_pstl.h

#ifndef AUTODIDACT_MODULES_ENGINE_SRC_PSTL_HPP_
#define AUTODIDACT_MODULES_ENGINE_SRC_PSTL_HPP_

#if defined(__GNUC__) && (__GNUC__ >= 9)

# include <execution>
# include <algorithm>
# include <numeric>
namespace exec = std::execution;

#elif defined(USE_LLVM_PSTL)

# include <__pstl_execution>
# include <__pstl_algorithm>
# include <__pstl_numeric>
namespace exec = std::execution;

#else

# include <pstl/execution>
# include <pstl/algorithm>
# include <pstl/numeric>

namespace exec = pstl::execution;

#endif

#endif //AUTODIDACT_MODULES_ENGINE_SRC_PSTL_HPP_
