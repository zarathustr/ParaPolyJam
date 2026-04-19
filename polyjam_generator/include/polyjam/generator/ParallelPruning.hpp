/*************************************************************************
 *                                                                       *
 * polyjam, a polynomial solver generator for C++                        *
 * Copyright (C) 2015 Laurent Kneip, The Australian National University  *
 *                                                                       *
 * Parallel row-basis pruning extension for Macaulay templates.           *
 *                                                                       *
 * This program is free software: you can redistribute it and/or modify  *
 * it under the terms of the GNU General Public License as published by  *
 * the Free Software Foundation, either version 3 of the License, or     *
 * (at your option) any later version.                                   *
 *                                                                       *
 *************************************************************************/

#ifndef POLYJAM_GENERATOR_PARALLELPRUNING_HPP_
#define POLYJAM_GENERATOR_PARALLELPRUNING_HPP_

#include <list>
#include <vector>

#include <polyjam/generator/CMatrix.hpp>

namespace polyjam
{
namespace generator
{
namespace pruning
{

std::vector<int> allRows( size_t rows );
std::list<int> toList( const std::vector<int> & rows );

// Selects a subset of the supplied original row indices whose rows form a
// row basis for that restricted matrix.  The selected rows are original rows,
// not linear-combination rows, so they can still be used as Polyjam template
// equations during code generation.
std::vector<int> rowBasisSequential(
    CMatrix & matrix,
    const std::vector<int> & rows );

// Divide-and-conquer row-basis pruning.  Each block is pruned independently,
// and block bases are merged in a reduction tree.  With oneTBB enabled the
// independent block reductions run in parallel; without oneTBB this function
// falls back to the same exact row-basis algorithm serially.
std::vector<int> rowBasisParallel(
    CMatrix & matrix,
    const std::vector<int> & rows,
    int requestedThreads,
    int minBlockRows );

int availableThreads();

}
}
}

#endif /* POLYJAM_GENERATOR_PARALLELPRUNING_HPP_ */
