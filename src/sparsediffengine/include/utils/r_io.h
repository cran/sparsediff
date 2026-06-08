/*
 * Copyright 2026 Daniel Cederberg and William Zhang
 *
 * This file is part of the SparseDiffEngine project.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
/*
 * r_io.h -- console-I/O redirection for the R `sparsediff` package build.
 *
 * GUARDED by SPARSEDIFF_R_PRINT (mirrors the SPARSEDIFF_USE_RBLAS hook in
 * cblas_wrapper.h):
 *
 *   - undefined (every upstream / standalone engine build): this header is
 *     empty, so behavior is byte-identical to upstream and the fork stays
 *     trivially mergeable.
 *
 *   - defined (the R `sparsediff` package, which compiles the engine through R's
 *     own toolchain and passes -DSPARSEDIFF_R_PRINT): every console write is
 *     routed to R's console via Rprintf / REprintf, and the libc symbols
 *     `printf` (_printf/_puts/_putchar), `stderr` (___stderrp) and `exit`
 *     (_exit) are never referenced -- so R CMD check's "checking compiled code"
 *     finds nothing to flag in the engine objects.
 *
 * The engine's I/O surface is small and was verified exhaustively: `printf`
 * (statistics in problem.c; no return-value use), `fprintf(stderr, ...)` (error
 * paths only -- the ONLY use of stderr, always as the first argument, no
 * vfprintf/fflush/fputs and no non-stderr fprintf), and `exit` (fatal paths).
 * There is no `abort`, no `stdout` token, and no `std::cout` (pure C).  So:
 *   printf  -> Rprintf
 *   stderr  -> a stream that forwards to REprintf (funopen/fopencookie; a
 *              sentinel + wrapped fprintf where neither exists, e.g. Windows)
 *   exit    -> Rf_error (noreturn longjmp to R's top level, matching exit's
 *              noreturn contract; better than killing the R process)
 *
 * Pure C: usable from every engine translation unit.  Delivered by an
 * #include from utils/tracked_alloc.h (which every flagged source includes),
 * so no non-portable -include compiler flag is needed.
 */
#ifndef SPARSEDIFFENGINE_UTILS_R_IO_H
#define SPARSEDIFFENGINE_UTILS_R_IO_H

#if defined(SPARSEDIFF_R_PRINT)

/* fopencookie (glibc) needs _GNU_SOURCE before any libc header. */
#if defined(__linux__) && !defined(_GNU_SOURCE)
#define _GNU_SOURCE 1
#endif

/* Pull in the real declarations BEFORE redefining printf / stderr / exit. */
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#if defined(__linux__)
#include <sys/types.h> /* ssize_t for the fopencookie write callback */
#endif

#include <R_ext/Error.h> /* Rf_error (noreturn) */
#include <R_ext/Print.h> /* Rprintf, REprintf, REvprintf */

/* ----- stderr -> a real FILE* whose writes go to REprintf --------------------
 * Using a real stream keeps fprintf(stderr, ...) correct without naming the
 * libc `stderr` symbol.  Declared static inline so it is shareable from every
 * TU and does not warn when a TU includes the header but uses no stderr. */
#if defined(__APPLE__) || defined(__FreeBSD__) || defined(__NetBSD__) || \
    defined(__OpenBSD__) || defined(__DragonFly__)

static int sparsediff_r_ewrite(void *cookie, const char *buf, int n) {
  (void)cookie;
  if (n > 0) REprintf("%.*s", n, buf);
  return n;
}
static inline FILE *sparsediff_r_stderr(void) {
  static FILE *f = NULL;
  if (!f) {
    f = funopen(NULL, NULL, &sparsediff_r_ewrite, NULL, NULL);
    if (f) setvbuf(f, NULL, _IONBF, 0);
  }
  return f;
}
#undef stderr
#define stderr (sparsediff_r_stderr())

#elif defined(__GLIBC__)

static ssize_t sparsediff_r_ewrite(void *cookie, const char *buf, size_t n) {
  (void)cookie;
  if (n > 0) REprintf("%.*s", (int)n, buf);
  return (ssize_t)n;
}
static inline FILE *sparsediff_r_stderr(void) {
  static FILE *f = NULL;
  if (!f) {
    cookie_io_functions_t io;
    memset(&io, 0, sizeof(io));
    io.write = &sparsediff_r_ewrite;
    f = fopencookie(NULL, "w", io);
    if (f) setvbuf(f, NULL, _IONBF, 0);
  }
  return f;
}
#undef stderr
#define stderr (sparsediff_r_stderr())

#else /* no funopen / fopencookie (e.g. Windows): sentinel + wrapped fprintf */

static char sparsediff_r_stderr_storage_;
#define SPARSEDIFF_R_STDERR ((FILE *)&sparsediff_r_stderr_storage_)
static inline int sparsediff_r_fprintf(FILE *stream, const char *fmt, ...) {
  va_list ap;
  int ret = 0;
  va_start(ap, fmt);
  if (stream == SPARSEDIFF_R_STDERR)
    REvprintf(fmt, ap);
  else
    ret = vfprintf(stream, fmt, ap);
  va_end(ap);
  return ret;
}
#undef stderr
#define stderr SPARSEDIFF_R_STDERR
#undef fprintf
#define fprintf(...) sparsediff_r_fprintf(__VA_ARGS__)

#endif

/* ----- printf -> Rprintf; exit -> Rf_error (done AFTER the system headers) -- */
#undef printf
#define printf(...) Rprintf(__VA_ARGS__)

#undef exit
#define exit(code) Rf_error("sparsediff engine: internal exit(%d)", (int)(code))

#endif /* SPARSEDIFF_R_PRINT */

#endif /* SPARSEDIFFENGINE_UTILS_R_IO_H */
