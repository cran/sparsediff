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
#ifndef CBLAS_WRAPPER_H
#define CBLAS_WRAPPER_H

#if defined(SPARSEDIFF_USE_RBLAS)
/* Portable build (e.g. the R 'sparsediff' package): route CBLAS calls through a
 * shim over the host's Fortran BLAS, so the engine uses the same BLAS as the
 * host application and stays portable to platforms whose BLAS provides no CBLAS
 * interface. Define -DSPARSEDIFF_USE_RBLAS and supply sparsediff_cblas.h. */
#include "sparsediff_cblas.h"
#elif defined(__APPLE__)
#define ACCELERATE_NEW_LAPACK
#include <Accelerate/Accelerate.h>
#else
#include <cblas.h>
#endif

#endif /* CBLAS_WRAPPER_H */
