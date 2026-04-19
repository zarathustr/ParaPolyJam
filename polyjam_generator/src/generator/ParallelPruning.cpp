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

#include <polyjam/generator/ParallelPruning.hpp>

#include <algorithm>
#include <cstdlib>
#include <iostream>

#ifndef POLYJAM_HAVE_TBB
#define POLYJAM_HAVE_TBB 0
#endif

#if POLYJAM_HAVE_TBB
#include <tbb/blocked_range.h>
#include <tbb/global_control.h>
#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>
#endif

#ifndef POLYJAM_PARALLEL_PRUNING_THREADS
#define POLYJAM_PARALLEL_PRUNING_THREADS 0
#endif

#ifndef POLYJAM_PARALLEL_PRUNING_MIN_BLOCK_ROWS
#define POLYJAM_PARALLEL_PRUNING_MIN_BLOCK_ROWS 64
#endif

namespace polyjam
{
namespace generator
{
namespace pruning
{

std::vector<int>
allRows( size_t rows )
{
  std::vector<int> result;
  result.reserve(rows);
  for( size_t i = 0; i < rows; ++i )
    result.push_back(static_cast<int>(i));
  return result;
}

std::list<int>
toList( const std::vector<int> & rows )
{
  std::list<int> result;
  for( size_t i = 0; i < rows.size(); ++i )
    result.push_back(rows[i]);
  return result;
}

int
availableThreads()
{
#if POLYJAM_HAVE_TBB
  int threads = static_cast<int>(tbb::task_scheduler_init::default_num_threads());
  if( threads > 0 )
    return threads;
#endif
  return 1;
}

static int
resolvedThreadCount( int requestedThreads )
{
  int threads = requestedThreads;
  if( threads <= 0 )
    threads = POLYJAM_PARALLEL_PRUNING_THREADS;
  if( threads <= 0 )
    threads = availableThreads();
  if( threads <= 0 )
    threads = 1;
#if !POLYJAM_HAVE_TBB
  threads = 1;
#endif
  return threads;
}

static int
resolvedMinBlockRows( int minBlockRows )
{
  if( minBlockRows <= 0 )
    minBlockRows = POLYJAM_PARALLEL_PRUNING_MIN_BLOCK_ROWS;
  if( minBlockRows <= 0 )
    minBlockRows = 64;
  return minBlockRows;
}

static void
deleteRows( std::vector<CMatrix::crow_t*> & rows )
{
  for( size_t i = 0; i < rows.size(); ++i )
    delete rows[i];
  rows.clear();
}

std::vector<int>
rowBasisSequential(
    CMatrix & matrix,
    const std::vector<int> & rows )
{
  std::vector<int> selectedRows;
  if( rows.empty() || matrix.cols() == 0 )
    return selectedRows;

  const int cols = static_cast<int>(matrix.cols());

  std::vector<CMatrix::crow_t*> work;
  std::vector<int> origins;
  work.reserve(rows.size());
  origins.reserve(rows.size());

  for( size_t i = 0; i < rows.size(); ++i )
  {
    CMatrix::crow_t * newRow = new CMatrix::crow_t();
    newRow->reserve(cols);
    bool nonzero = false;
    for( int c = 0; c < cols; ++c )
    {
      core::Coefficient coeff = matrix(rows[i], c);
      if( !coeff.isZero() )
        nonzero = true;
      newRow->push_back(coeff);
    }

    if( nonzero )
    {
      work.push_back(newRow);
      origins.push_back(rows[i]);
    }
    else
    {
      delete newRow;
    }
  }

  if( work.empty() )
    return selectedRows;

  core::Coefficient zero = (*(work[0]))[0].zero();
  core::Coefficient one = (*(work[0]))[0].one();

  const int rowCount = static_cast<int>(work.size());
  int frontRow = 0;
  int currentCol = 0;

  while( frontRow < rowCount && currentCol < cols )
  {
    int pivotRow = -1;
    for( int r = frontRow; r < rowCount; ++r )
    {
      if( (*(work[r]))[currentCol] != zero )
      {
        pivotRow = r;
        break;
      }
    }

    if( pivotRow < 0 )
    {
      ++currentCol;
      continue;
    }

    if( pivotRow != frontRow )
    {
      std::swap(work[pivotRow], work[frontRow]);
      std::swap(origins[pivotRow], origins[frontRow]);
    }

    selectedRows.push_back(origins[frontRow]);

    core::Coefficient leadingCoefficient = (*(work[frontRow]))[currentCol] + zero;
    (*(work[frontRow]))[currentCol] = one + zero;
    core::Coefficient leadingCoefficientInv = (*(work[frontRow]))[currentCol] / leadingCoefficient;

    for( int c = currentCol + 1; c < cols; ++c )
      (*(work[frontRow]))[c] *= leadingCoefficientInv;

    std::vector<int> nonzeroColumns;
    nonzeroColumns.reserve(cols - currentCol);
    for( int c = currentCol; c < cols; ++c )
    {
      if( (*(work[frontRow]))[c] != zero )
        nonzeroColumns.push_back(c);
    }

    for( int r = frontRow + 1; r < rowCount; ++r )
    {
      core::Coefficient multiplier = (*(work[r]))[currentCol] + zero;
      if( multiplier != zero )
      {
        for( size_t ci = 0; ci < nonzeroColumns.size(); ++ci )
        {
          const int c = nonzeroColumns[ci];
          (*(work[r]))[c] -= multiplier * (*(work[frontRow]))[c];
        }
      }
    }

    ++frontRow;
    ++currentCol;
  }

  deleteRows(work);
  return selectedRows;
}

static std::vector< std::vector<int> >
partitionRows( const std::vector<int> & rows, int blockCount )
{
  std::vector< std::vector<int> > blocks;
  if( rows.empty() )
    return blocks;

  if( blockCount < 1 )
    blockCount = 1;
  if( blockCount > static_cast<int>(rows.size()) )
    blockCount = static_cast<int>(rows.size());

  blocks.resize(blockCount);
  const size_t n = rows.size();
  for( int b = 0; b < blockCount; ++b )
  {
    const size_t begin = (n * static_cast<size_t>(b)) / static_cast<size_t>(blockCount);
    const size_t end   = (n * static_cast<size_t>(b + 1)) / static_cast<size_t>(blockCount);
    blocks[b].reserve(end - begin);
    for( size_t i = begin; i < end; ++i )
      blocks[b].push_back(rows[i]);
  }
  return blocks;
}

static std::vector<int>
concatRows( const std::vector<int> & a, const std::vector<int> & b )
{
  std::vector<int> merged;
  merged.reserve(a.size() + b.size());
  merged.insert(merged.end(), a.begin(), a.end());
  merged.insert(merged.end(), b.begin(), b.end());
  return merged;
}

std::vector<int>
rowBasisParallel(
    CMatrix & matrix,
    const std::vector<int> & rows,
    int requestedThreads,
    int minBlockRows )
{
  if( rows.empty() )
    return rows;

  const int threads = resolvedThreadCount(requestedThreads);
  const int minRows = resolvedMinBlockRows(minBlockRows);

  if( threads <= 1 || static_cast<int>(rows.size()) <= minRows )
    return rowBasisSequential(matrix, rows);

  int blockCount = threads * 4;
  if( blockCount < 1 )
    blockCount = 1;
  const int maxUsefulBlocks = std::max(1, static_cast<int>((rows.size() + minRows - 1) / minRows));
  blockCount = std::min(blockCount, maxUsefulBlocks);
  blockCount = std::min(blockCount, static_cast<int>(rows.size()));

  if( blockCount <= 1 )
    return rowBasisSequential(matrix, rows);

  std::vector< std::vector<int> > blocks = partitionRows(rows, blockCount);

#if POLYJAM_HAVE_TBB
  tbb::global_control control(tbb::global_control::max_allowed_parallelism, threads);
  tbb::parallel_for(
      tbb::blocked_range<int>(0, static_cast<int>(blocks.size())),
      [&]( const tbb::blocked_range<int> & range )
      {
        for( int b = range.begin(); b != range.end(); ++b )
          blocks[b] = rowBasisSequential(matrix, blocks[b]);
      });
#else
  for( int b = 0; b < static_cast<int>(blocks.size()); ++b )
    blocks[b] = rowBasisSequential(matrix, blocks[b]);
#endif

  while( blocks.size() > 1 )
  {
    const int pairCount = static_cast<int>(blocks.size() / 2);
    const bool hasOdd = (blocks.size() % 2) != 0;
    std::vector< std::vector<int> > next(pairCount + (hasOdd ? 1 : 0));

#if POLYJAM_HAVE_TBB
    tbb::parallel_for(
        tbb::blocked_range<int>(0, pairCount),
        [&]( const tbb::blocked_range<int> & range )
        {
          for( int p = range.begin(); p != range.end(); ++p )
          {
            std::vector<int> merged = concatRows(blocks[2 * p], blocks[2 * p + 1]);
            next[p] = rowBasisSequential(matrix, merged);
          }
        });
#else
    for( int p = 0; p < pairCount; ++p )
    {
      std::vector<int> merged = concatRows(blocks[2 * p], blocks[2 * p + 1]);
      next[p] = rowBasisSequential(matrix, merged);
    }
#endif

    if( hasOdd )
      next.back() = blocks.back();

    blocks.swap(next);
  }

  if( blocks.empty() )
    return std::vector<int>();
  return blocks.front();
}

}
}
}
