/*****************************************************************************
 *   GATB : Genome Assembly Tool Box
 *   Copyright (C) 2014  INRIA
 *   Authors: G.Rizk, G.Benoit
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*****************************************************************************/

#ifndef _BLOOCOO_HPP_
#define _BLOOCOO_HPP_

/********************************************************************************/

#include <gatb/gatb_core.hpp>

#include <string>
#include <sstream>

#define BLOOCOO_VERSION_MAJOR 1
#define BLOOCOO_VERSION_MINOR 0
#define BLOOCOO_VERSION_PATCH 6


#include "ReadWriter.h"
//#include <vector>
//#include <bitset>

using namespace std;
#define TAB_MULTIVOTE_SIZE (16*16*128)

/********************************************************************************/
/* Class Bloocoo for read correction : takes as input a bag of solid kmers (dsk result),
 insert that in a bloom filter, and correct read form it*/
/********************************************************************************/

typedef Kmer<>::ModelDirect    ModelDirect;
typedef Kmer<>::ModelCanonical ModelCanonical;
typedef Kmer<>::Type  kmer_type;
typedef Kmer<>::Count kmer_count;
//typedef gatb::core::tools::collections::impl OaHash;


class Bloocoo : public Tool
{
//private:
public:

	enum CorrectionMethod{
		TWO_SIDED,
		AGRESSIVE,
		VOTE,
		MULTI_MUTATE_VOTE,
		NB_CORRECTION_METHODS,
	};
	
	/** Synchronizer. */
	ISynchronizer* errtab_mutex;
	ISynchronizer* getErrtabMutex()
	{
	    if (errtab_mutex==0)  { errtab_mutex = System::thread().newSynchronizer(); }
	    return errtab_mutex;
	}
    
    std::string _solidFile;
   // uint64_t        _seq_num; // removed, was not thread safe
    size_t _kmerSize;
    int _nb_kmers_checked;
    int _nb_min_valid;
	int _nb_less_restrictive_correction;
    int _max_multimutation_distance;
    int _only_decrease_nb_min_valid;
    
    IBank* _inputBank;

    static const char* STR_ION;
    //static const char* STR_NB_MIN_VALID;
    //static const char* STR_NB_VALIDATED_KMERS;
    static const char* STR_ERR_TAB;
	static const char* STR_MANUAL_SOLIDITY;
	static const char* STR_RECALL;
	static const char* STR_HISTO_ONLY;
	static const char* STR_FROM_H5;
	static const char* STR_PRECISION;
	static const char* STR_SLOW;
	static const char* STR_MAX_TRIM;
	static const char* STR_BIT_BLOOM_PER_KMER;

	
    bool _wantErrTabFile;
    IFile*      _errfile;
    IFile*      _errfile_full;
    IFile*      _debug;

	
	bool _fromH5Mode ;
	bool _countOnlyMode ;
	bool _manualSolidity;
	
	bool _ion_mode;
    unsigned int _max_trim;
    
	int __correction_methods_successes[NB_CORRECTION_METHODS]; //Variable de debug pour afficher l'efficacité des méthodes de correction
	int __correction_methods_calls[NB_CORRECTION_METHODS]; //Variable de debug pour afficher le nb d'appels
	

	
public:

    Bloocoo ();
    ~Bloocoo ();

	
private:

    IBloom<kmer_type>* _bloom;
    
    void execute ();
    virtual IBloom<kmer_type>* createBloom ();
    void chooseCorrectionParams();
    
};


/********************************************************************************/
/* Class CorrectReads used to correct a single read. Bloocoo creates CorrectReads instances
 * to perform the correction in multiple threads. Each instances of this classes
 * is synchronized*/
/********************************************************************************/
class CorrectReads
{
	public:
	
		CorrectReads(IBloom<kmer_type>* bloom, Bag<Sequence>* outbank, Bloocoo * bloocoo, u_int64_t* nb_errors_corrected, u_int64_t* nb_ins_corrected, u_int64_t* nb_del_corrected, int nb_cores, int * nbliving);
		CorrectReads(const CorrectReads& cr);
		~CorrectReads ();
		void operator()(Sequence& sequence);
		
		enum Direction //typedef enum direction
		{
			LEFT,
			RIGHT
		};// direction_t;
		
	private:
	
		IBloom<kmer_type> * _bloom; // the bloom containing the solid kmers
		Bag<Sequence> *             _outbank; // the bank cto insert the result : corrected reads
		Bloocoo *         _bloocoo; // the parent bloocoo object
		unsigned char *   _tab_multivote;
		int *  _nb_living;
		int _thread_id;
		char * _tempseq;
		bool _use_newSeq;
		Sequence * _newSeq; //used for indel correction

		OrderedBankWriter * _bankwriter;
		void setBankWriter (OrderedBankWriter* bankwriter) { SP_SETATTR(bankwriter); }
		
		Sequence* _sequence;
		char* _readseq;
		vector<int> _corrected_pos;
		vector<kmer_type> _kmers;
		
		size_t _kmerSize;
		int _readSize;
		int _kmerCount;
		int _nb_kmers_checked;
		int _nb_min_valid;
		int _nb_kmer_offset;
		int _nb_less_restrictive_correction;
		int _max_multimutation_distance;
		bool _only_decrease_nb_min_valid;
		bool _wantErrTabFile;
		
		bool _continue_correction;
		
		// KmerModel           model;
		// KmerModel::Iterator itKmer;

		ISynchronizer* _synchro;
		void setSynchro (ISynchronizer* synchro)  { SP_SETATTR(synchro); }
		ISynchronizer* getSynchro ()  { return _synchro; }

		u_int64_t *  _total_nb_errors_corrected;
		u_int64_t *  _total_nb_ins_corrected;
		u_int64_t *  _total_nb_del_corrected;

		u_int64_t   _local_nb_errors_corrected;
		u_int64_t   _local_nb_ins_corrected;
		u_int64_t   _local_nb_del_corrected;

		size_t _temp_nb_seq_done;
		
		void execute();
		void writeSequence();
		void trimSequence();
		void update_nb_errors_corrected(int nb_errors_corrected); //not needed if execute2 is used
		int apply_correction(int pos, int good_nt,int algo);
		int twoSidedCorrection(int pos);
		int aggressiveCorrection(int pos, int nb_kmer_check, Direction direction);
		int voteCorrectionInUntrustedZone(int start_pos, int end_pos, int nb_kmer_checked);
		int voteCorrection(int start_pos, int nb_kmer_check);

		int multiMutateVoteCorrectionInUntrustedZone(int start_pos, int end_pos, int nb_kmer_checked,int expected_first_pos,int expected_second_pos);
		int multiMutateVoteCorrection(int start_pos, int nb_kmer_check,int expected_first_pos,int expected_second_pos);
		void multiMutateVoteCorrectionRec(int start_pos, int kmer_offset, kmer_type current_kmer, int nb_kmer_check, int kmer_index, int current_nb_mutation, int* max_vote, int* nb_max_vote, int *good_index, int pos1, int nt1,int expected_first_pos,int expected_second_pos);
		
		void extendedAgressiveCorrection(int votes[4], ModelCanonical* model, kmer_type* last_mutated_kmer, int mutated_nt, int nb_kmer_check, Direction direction);
		void extendedAgressiveCorrectionRec(int votes[4], ModelCanonical* model, kmer_type* mutated_kmer, int mutated_nt, int nb_kmer_check, Direction direction, bool posVoted[], int depth);
		
		void codeSeedBin(ModelCanonical* model, kmer_type* kmer, int nt, Direction direction);
		void codeSeedNT(ModelCanonical* model, kmer_type* kmer, char nt, Direction direction);
		
		//void print_agressive_votes(int votes[4]);
		//void print_votes(int votes[][4], int nb_column);
		//void print_read_correction_state(KmerModel* model);
		//void print_read_if_not_fully_corrected(KmerModel* model);
		void mutate_kmer(kmer_type * kmer, int pos, char nt);
		bool is_pos_correctable(int pos);
	
		unsigned int make_index(int pos1, int dist, int nt1, int nt2);
		int decode_index(unsigned int idx, int * pos1, int * dist, int * nt1, int * nt2);
		
		
		
		
		
		//ion correction
		//kmer_type _last_wrong_kmer;
		//int _pos_homopo;
		//void ionCorrection(int ii);
		//void ionCorrection2(int ii, kmer_type current_kmer, int untrusted_zone_size);
		
		
		
		//new method
		void execute2();
		bool searchError(bool* error_exist);
		int searchErrorZoneRec(int pos, bool extendRight, int trustedKmerCount);
		int startCorrectionInZone(int startPos, int endPos);
		
		void executeIonCorrection();
		

		
};

/********************************************************************************/

#endif /* _BLOOCOO_HPP_ */

