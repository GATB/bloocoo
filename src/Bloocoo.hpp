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
#include <sstream>

#define TAB_MULTIVOTE_SIZE (16*16*128)

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

	enum Direction //typedef enum direction
	{
		LEFT,
		RIGHT
	};// direction_t;

	enum CorrectionMethod{
		TWO_SIDED,
		AGRESSIVE,
		VOTE,
		MULTI_MUTATE_VOTE,
		READ_SIDE,
		NB_CORRECTION_METHODS,
	};
	
	
    size_t          _kmerSize;
    std::string     _solidFile;
   // uint64_t        _seq_num; // removed, was not thread safe
    int            _nb_kmers_checked;
    
    int            _nb_min_valid;
    
    IBank* _inputBank;

    static const char* STR_NB_MIN_VALID;
    static const char* STR_NB_VALIDATED_KMERS;
    
    IFile*      _errfile;
    IFile*      _errfile_full;
    IFile*      _debug;

	
	std::string __badReadStack; //Variable de debug qui contient l'empreinte de la correction des reads
	int __correction_methods_successes[NB_CORRECTION_METHODS]; //Variable de debug pour afficher l'efficacité des méthodes de correction
	int __correction_methods_calls[NB_CORRECTION_METHODS]; //Variable de debug pour afficher le nb d'appels
	
	int _max_multimutation_distance;
//	std::ostringstream _oss; //Use to convert int to string many time per frame
//	std::istringstream _iss; //Use to convert string to int many time per frame
	
public:

    /** */
    Bloocoo ();

	
    unsigned int make_index(int pos1, int dist, int nt1, int nt2);
    int decode_index(unsigned int idx, int * pos1, int * dist, int * nt1, int * nt2);

    
	void update_nb_errors_corrected(int nb_errors_corrected, u_int64_t* _local_nb_errors_corrected, bool* continue_correction);
	int apply_correction(char *readseq, int pos, int good_nt,int algo, Sequence& cur_seq,std::vector<int>& corrected_pos);
	int twoSidedCorrection(int pos, char *readseq, kmer_type* kmers[], Sequence& s,std::vector<int>& corrected_pos);
	int aggressiveCorrection(int pos, char *readseq, kmer_type* kmers[], int nb_kmer_check, int readlen, Direction direction, Sequence& s,std::vector<int>& corrected_pos, int min_nb_kmer_valid_offset);
	int voteCorrectionInUntrustedZone(int start_pos, int end_pos, char *readseq, kmer_type* kmers[], int nb_kmer_checked, Sequence& s,std::vector<int>& corrected_pos, int min_nb_kmer_valid_offset);
	int voteCorrection(int start_pos, char *readseq, kmer_type* kmers[], int nb_kmer_check, Sequence& s,std::vector<int>& corrected_pos, int min_nb_kmer_valid_offset);

	int multiMutateVoteCorrectionInUntrustedZone(int start_pos, int end_pos, char *readseq, kmer_type* kmers[], int nb_kmer_checked, unsigned char* _tab_multivote,int expected_first_pos,int expected_second_pos, int readlen, Sequence& s,std::vector<int>& corrected_pos, int min_nb_kmer_valid_offset);
	int multiMutateVoteCorrection(int start_pos, char *readseq, kmer_type* kmers[], int nb_kmer_check, unsigned char* _tab_multivote,int expected_first_pos,int expected_second_pos, int readlen, Sequence& s,std::vector<int>& corrected_pos, int min_nb_kmer_valid_offset);
	int multiMutateVoteCorrectionRec(int start_pos, int kmer_offset, kmer_type current_kmer, char *readseq, kmer_type* kmers[], int nb_kmer_check, int kmer_index, int current_nb_mutation, unsigned char* _tab_multivote, int* max_vote, int* nb_max_vote, int *good_index, int pos1, int nt1,int expected_first_pos,int expected_second_pos, Sequence& s,std::vector<int>& corrected_pos);
	//int multiMutateVoteCorrectionRec(int start_pos, int kmer_offset, kmer_type current_kmer, char *readseq, kmer_type* kmers[], int nb_kmer_check, int kmer_index, int current_nb_mutation, unsigned char* _tab_multivote, int pos1, int nt1);
	
	void extendedAgressiveCorrection(int votes[4], KmerModel* model, kmer_type* last_mutated_kmer, int mutated_nt, int nb_kmer_check, Direction direction);
	void extendedAgressiveCorrectionRec(int votes[4], KmerModel* model, kmer_type* mutated_kmer, int mutated_nt, int nb_kmer_check, Direction direction, bool posVoted[], int depth);
	
	void codeSeedBin(KmerModel* model, kmer_type* kmer, int nt, Direction direction);
	void codeSeedNT(KmerModel* model, kmer_type* kmer, char nt, Direction direction);
	
	void print_agressive_votes(int votes[4]);
	void print_votes(int votes[][4], int nb_column);
	void print_read_correction_state(KmerModel* model, Sequence& s);
	void print_read_if_not_fully_corrected(KmerModel* model, Sequence& s);
	void mutate_kmer(kmer_type * kmer, int pos, char nt);
	bool is_pos_correctable(int pos, char* readseq,std::vector<int>& corrected_pos);

	
private:

    collections::impl::Bloom<kmer_type>* _bloom;
    
    /** */
    void execute ();
    
    

    /** */
    virtual collections::impl::Bloom<kmer_type>* createBloom ();
    
    
    
};

/********************************************************************************/

#endif /* _BLOOCOO_HPP_ */

