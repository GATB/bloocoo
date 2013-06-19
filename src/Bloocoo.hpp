/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

#ifndef _BLOOCOO_HPP_
#define _BLOOCOO_HPP_

/********************************************************************************/

#include <gatb/tools/misc/impl/Tool.hpp>

#include <gatb/tools/math/Integer.hpp>
#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>
#include <gatb/tools/designpattern/impl/Command.hpp>
#include <gatb/tools/misc/impl/TimeInfo.hpp>
#include <gatb/tools/misc/impl/OptionsParser.hpp>
#include <gatb/bank/impl/Bank.hpp>
#include <gatb/bank/impl/BankBinary.hpp>
#include <gatb/kmer/impl/Model.hpp>
#include <gatb/tools/collections/impl/Bloom.hpp>

#include <string>

/********************************************************************************/

/** NOTE: we should not include namespaces here => only to make user life easier... */
using namespace gatb::core;
using namespace gatb::core::tools;
using namespace gatb::core::bank;
using namespace gatb::core::kmer::impl;

using namespace gatb::core::system;
using namespace gatb::core::system::impl;

/********************************************************************************/
/* Class Bloocoo for read correction : takes as input a bag of solid kmers (dsk result),
 insert that in a bloom filter, and correct read form it*/
/********************************************************************************/

class Bloocoo : public misc::impl::Tool
{
//private:
public:

    size_t          _kmerSize;
    std::string     _solidFile;
    uint64_t        _seq_num;
    int            _nb_kmers_checked;
    
    int            _nb_passes_per_read;
    
    IBank* _inputBank;

    static const char* STR_NB_ITER_PER_READ;
    static const char* STR_NB_VALIDATED_KMERS;
    
    IFile*      _errfile;

public:

    /** */
    Bloocoo ();

private:

    
    /** */
    void execute ();

    /** */
    virtual collections::impl::Bloom<kmer_type>* createBloom ();
};

/********************************************************************************/

#endif /* _BLOOCOO_HPP_ */

