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
	std::vector<int> _corrected_pos;
	
	std::string __badReadStack; //Variable de debug qui contient l'empreinte de la correction des reads
	int __correction_methods_successes[3]; //Variable de debug pour afficher l'efficacité des méthodes de correction
	
public:

    /** */
    Bloocoo ();

	
	enum Direction //typedef enum direction
	{
		LEFT,
		RIGHT
	};// direction_t;

	void update_nb_errors_corrected(int nb_errors_corrected, u_int64_t* _local_nb_errors_corrected, bool* continue_correction);
	int apply_correction(char *readseq, int pos, int good_nt);
	int twoSidedCorrection(int pos, char *readseq, kmer_type* kmers[]);
	int aggressiveCorrection(int pos, char *readseq, kmer_type* kmers[], int nb_kmer_check, Direction direction);
	int voteCorrectionInUntrustedZone(int start_pos, int end_pos, char *readseq, kmer_type* kmers[], int nb_kmer_checked);
	int voteCorrection(int start_pos, char *readseq, kmer_type* kmers[], int nb_kmer_check);

	kmer_type codeSeedBin(KmerModel* model, kmer_type* kmer, int nt, Direction direction);
	kmer_type codeSeedNT(KmerModel* model, kmer_type* kmer, char nt, Direction direction);
	
	void print_agressive_votes(int votes[4]);
	void print_votes(int votes[][4], int nb_column);
	void print_read_correction_state(KmerModel* model, Sequence& s, int seq_num);
	void print_read_if_not_fully_corrected(KmerModel* model, Sequence& s);
	void mutate_kmer(kmer_type * kmer, int pos, char nt);
	bool is_pos_correctable(int pos, char* readseq);

	
private:

    collections::impl::Bloom<kmer_type>* _bloom;
    
    /** */
    void execute ();
    
    

    /** */
    virtual collections::impl::Bloom<kmer_type>* createBloom ();
    
    
    
};

/********************************************************************************/

#endif /* _BLOOCOO_HPP_ */

