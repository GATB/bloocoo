/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

#include <Bloocoo.hpp>

#include <libgen.h>

#include "ReadWriter.h"

// We use the required packages
using namespace std;

#define DEBUG(a)  printf a

/********************************************************************************/
// We define some string constants.
/********************************************************************************/
const char* Bloocoo::STR_NB_MIN_VALID         = "-nbmin-valid";
const char* Bloocoo::STR_NB_VALIDATED_KMERS         = "-nkmer-checked";


// these 2 tabs are now known globally
//char bin2NT[4] = {'A','C','T','G'};
//char binrev[4]    = {2,3,0,1};

//char bin2NTrev[4] = {'T','G','A','C'};



/********************************************************************************/
//fonctions for correction
/********************************************************************************/



//PRINT_DEBUG = 1: print all reads states
//PRINT_DEBUG = 2: print reads states not fully corrected
#define PRINT_DEBUG 0

#define PRINT_STATS 0 // only in single thread mode

#define PRINT_LOG_MULTI 0


//#define SERIAL // to force serial mode









/********************************************************************************/
//functor for read correction this is the code to correct a single read
/********************************************************************************/
class CorrectReads
{
public:
    void operator() ( Sequence& current_seq) //no const otherwise error with tKmer.setData
    {

//    	_bloocoo->_seq_num = s.getIndex();
    	//printf("SEQ NUM: %i\n", current_seq.getIndex());
        //if(s.getIndex() != 681){
        //	return;
        //}
        
        // _bloom   is the bloom filter containing the solid kmers
        // _bloom->contains(kmer)  to ask if it contains a kmer
        
        _corrected_pos.clear();
        
        if(PRINT_DEBUG){ _bloocoo->__badReadStack = "\n\n\n"; }
        

        
        //printf("---- new read ----\n");
        
        char * readseq = current_seq.getDataBuffer(); // the nucleotide sequence of the read
        size_t sizeKmer = _bloocoo->_kmerSize;
        
        //int nbPasses = _bloocoo->_nb_passes_per_read;
        int max_nb_kmers_checked = _bloocoo->_nb_kmers_checked;
        int readlen = current_seq.getDataSize();
        
        int nb_checked;
        KmerModel model (sizeKmer);
        KmerModel::Iterator itKmer (model);
		kmer_type current_kmer;
		kmer_type current_kmer_min;
			
		int nb_kmer_offset = -2;
        
        for(int j=0; j<_nb_less_restrictive_correction; j++){
			nb_kmer_offset += 2;
			
			//printf("%i\n", nb_kmer_offset);
			
			int nb_kmer_trusted = 0;
			bool continue_correction = true;
			bool first_gap = true;
			int ii=0;
			
			
			while(continue_correction)
			{

				nb_kmer_trusted = 0;
				continue_correction = false;
				first_gap = true;
				ii = 0;
				
				itKmer.setData (current_seq.getData(),KMER_DIRECT);
				
				int untrusted_zone_size = 0;
				int trusted_zone_size = 0;
				int real_untrusted_zone_size = 0;
				
				//Mettre en dehors du while dans une version final (attention dangereux), faire gaffe a ne jamais utiliser un kmer de ce tableau
				//avec un indice > ii
				kmer_type* kmers[readlen-sizeKmer+1]; 
														
							
				if(PRINT_DEBUG){ _bloocoo->print_read_correction_state(&model, current_seq); }

			   // getSynchro()->lock()  ;

				// We iterate the kmers of this sequence
				for (itKmer.first(); !itKmer.isDone(); itKmer.next(),ii++)
				{
					
				   

					int nb_errors_cor = 0;
					
					current_kmer = *itKmer;
					kmers[ii] = &(*itKmer);
					
					current_kmer_min = min(revcomp(current_kmer, sizeKmer), current_kmer);

					//Si on veut checker le dernier kmer du read et qu'on est dans une zone untrusted alors on considere
					//automatiquement ce dernier kmer comme untrusted et on laisse les démarches de correction de fin de read s'effectuer.
					bool is_last_kmer_indexed_after_hole = (ii==readlen-sizeKmer && trusted_zone_size==0);

					
					if (!is_last_kmer_indexed_after_hole && _bloom->contains(current_kmer_min)) //kmer is solid
					{
						
						nb_kmer_trusted += 1;
						
						trusted_zone_size += 1;
						
						//beginning of indexed zone
						if(trusted_zone_size == 2) 
						{
						   
							if (untrusted_zone_size>1){
							
								//if(ii-2 == 30){
								//	printf("%i %i %i %i %i\n", (untrusted_zone_size-1) == sizeKmer, ii > sizeKmer, ii-2, untrusted_zone_size, sizeKmer);
								//}
								
								//two-sided conservative correction
								// this should be an isolated error, middle of the read
								// if error at pos sizeKmercorrected_pos+1+ii, could work in theory but would need kmer_begin+1
								if(((untrusted_zone_size-1) == sizeKmer)  &&  (ii > sizeKmer+2)){
									
									//if(ii-2==29){
									//	printf("%i %i\n", untrusted_zone_size, ii-2);
									//}
									
									if(PRINT_DEBUG){ _bloocoo->__badReadStack += "\t\tTwo sided (hole size k)\n";}
									nb_errors_cor = _bloocoo->twoSidedCorrection(ii-2, readseq, kmers,current_seq,_corrected_pos);

								}
							
								if(nb_errors_cor == 0){
									nb_checked = max_nb_kmers_checked;
									nb_checked = min(max_nb_kmers_checked, untrusted_zone_size-1); // -2 changed to -1 : was caus of pb for sides
									nb_checked = max(nb_checked, 0);
									
									//if first gap, cannot correct from the left side of the gap
									//2eme condition: Agressive right peut corriger le debut si le trou est plus long que sizeKmer mais est-ce utile?
									if(!first_gap)// || (first_gap && untrusted_zone_size > sizeKmer))
									{
										if(PRINT_DEBUG){ _bloocoo->__badReadStack += "\t\tAggressive correction (big hole 1 right)\n"; }
										nb_errors_cor = _bloocoo->aggressiveCorrection(ii-untrusted_zone_size, readseq,  kmers, nb_checked, readlen, Bloocoo::RIGHT,current_seq,_corrected_pos, nb_kmer_offset);
									}
									
									
									if(nb_errors_cor == 0){
										if(PRINT_DEBUG){ _bloocoo->__badReadStack += "\t\tAggressive correction (big hole 2 left)\n"; }
										nb_errors_cor = _bloocoo->aggressiveCorrection(ii-2 , readseq,  kmers, nb_checked, readlen, Bloocoo::LEFT,current_seq,_corrected_pos, nb_kmer_offset);
									}
									
									if(nb_errors_cor == 0){
										if(PRINT_DEBUG){ _bloocoo->__badReadStack += "\t\tVote correction (big hole 2)\n"; }
										nb_errors_cor = _bloocoo->voteCorrectionInUntrustedZone(ii- untrusted_zone_size, ii-2, readseq, kmers, nb_checked,current_seq,_corrected_pos, nb_kmer_offset);
									}

									if(nb_errors_cor == 0){
										if(PRINT_DEBUG){ _bloocoo->__badReadStack += "\t\tMulti Mutate Vote correction (big hole 2)\n"; }
										nb_errors_cor = _bloocoo->multiMutateVoteCorrectionInUntrustedZone(ii- untrusted_zone_size, ii-2, readseq, kmers, nb_checked, _tab_multivote,-1,ii-2,readlen,current_seq,_corrected_pos, nb_kmer_offset);
									}
									
								
								}
							}
							
						}
						
						if(trusted_zone_size==2){
							first_gap = false;
						}
						
						////////////traitement d'un faux positif du bloom isolé
						//   if(trusted_zone_size==1) tai_previous_break = untrusted_zone_size; // begin of indexed zone    ( if untrusted_zone_size ?)
						if (trusted_zone_size > 1) // do not reset  untrusted_zone_size if a single positive (might be a FP)
						{
							untrusted_zone_size = 0; //reset previous hole size
						}
						if (trusted_zone_size==1)// begin of indexed zone
						{
							untrusted_zone_size++; //  treat a single solidkmer  as an erroneous kmer
						}
						//////////
						
						//version sans traitement du FP du bloom
						//if(trusted_zone_size==1){
						//	untrusted_zone_size = 0;
						//}
						
					}
					else //kmer contains an error
					{

						untrusted_zone_size += 1;
						trusted_zone_size = 0;
								
						//end of the read, we should treat this gap here correc snp with trad method
						if(ii == (readlen-sizeKmer)){
							
							/*
							//If the last kmer is indexed, we 
							if (_bloom->contains(current_kmer_min)){
								nb_kmer_trusted += 1;
								if(untrusted_zone_size-1==sizeKmer){
									if(PRINT_DEBUG){ _bloocoo->__badReadStack += "\t\tTwo sided (hole size k)\n";}
									nb_errors_cor = _bloocoo->twoSidedCorrection(ii-1, readseq, kmers,current_seq,_corrected_pos);
								}
							}
							//Apply twosided if the untrusted zone size is equal to kmerSize
							else if(untrusted_zone_size==sizeKmer){
								if(PRINT_DEBUG){ _bloocoo->__badReadStack += "\t\tTwo sided (hole size k)\n";}
								nb_errors_cor = _bloocoo->twoSidedCorrection(ii, readseq, kmers,current_seq,_corrected_pos);
								//printf("%s", _bloocoo->__badReadStack.c_str());
								//printf("%i\n", ii);
							}*/
							
							if(nb_errors_cor == 0){
								nb_checked = max_nb_kmers_checked;
								nb_checked =  min(max_nb_kmers_checked, untrusted_zone_size);
								nb_checked = max(nb_checked, 0);
								
								if(PRINT_DEBUG){ _bloocoo->__badReadStack += "\t\tAggressive correction (end)\n"; }
								nb_errors_cor = _bloocoo->aggressiveCorrection(readlen - untrusted_zone_size - sizeKmer + 1, readseq,  kmers, nb_checked, readlen ,Bloocoo::RIGHT,current_seq,_corrected_pos, nb_kmer_offset);
								
								if(nb_errors_cor == 0){
									if(PRINT_DEBUG){ _bloocoo->__badReadStack += "\t\tVote correction (end)\n"; }
									nb_errors_cor = _bloocoo->voteCorrectionInUntrustedZone(readlen - (untrusted_zone_size-1) - sizeKmer , readlen-sizeKmer, readseq, kmers, nb_checked,current_seq,_corrected_pos, nb_kmer_offset);
								}
								
								if(nb_errors_cor == 0){
									if(PRINT_DEBUG){ _bloocoo->__badReadStack += "\t\tMulti Mutate Vote correction (end)\n"; }
									nb_errors_cor = _bloocoo->multiMutateVoteCorrectionInUntrustedZone(readlen - (untrusted_zone_size-1) - sizeKmer, readlen-sizeKmer, readseq, kmers, nb_checked, _tab_multivote,readlen - untrusted_zone_size ,-1,readlen,current_seq,_corrected_pos, nb_kmer_offset);
								}
							}
						}
						
					}
					
					_bloocoo->update_nb_errors_corrected(nb_errors_cor, &_local_nb_errors_corrected, &continue_correction);

				} // end of kmers iteration over the read
				

			}
			
			if(nb_kmer_trusted == readlen-sizeKmer+1){
				break;
			}
        }
        
        if(PRINT_DEBUG == 1){ 
        	printf("%s", _bloocoo->__badReadStack.c_str());
        }
        else if(PRINT_DEBUG == 2){
        	_bloocoo->print_read_if_not_fully_corrected(&model, current_seq);
        }
        

#ifdef SERIAL
        _outbank->insert (current_seq); //output corrected sequence
#else
        _bankwriter->insert (current_seq);
        _temp_nb_seq_done ++;
        
        if(_temp_nb_seq_done == 10000) // careful with this val, should be a divisor of buffer size in ordered bank
        {
            _bankwriter->incDone(_temp_nb_seq_done);
            _temp_nb_seq_done=0;
        }
#endif

       
        
        //    printf("%s\n",s.getDataBuffer());
        //    printf("%s\n",readseq);
        
       // _bloocoo->_seq_num++; //counter 0 based, not thread-safe  //todo : include the sequence number in the sequence type
        
    }

    
    //default constructor
    CorrectReads ()
    : _bloom(NULL), _outbank(NULL), _bloocoo(NULL),
    _total_nb_errors_corrected (NULL), _local_nb_errors_corrected(0),_bankwriter(NULL),_temp_nb_seq_done(0)
    {
       // printf("------- CorrectReads default Constructor  %p --------- \n",this);
        _tab_multivote = (unsigned char *) malloc(TAB_MULTIVOTE_SIZE*sizeof(unsigned char)); // pair of muta  = 16 nt *128 pos * 16 (max dist)
    }
    
    
    CorrectReads (Bloom<kmer_type>* bloom, Bag<Sequence>* outbank, Bloocoo * bloocoo, u_int64_t* nb_errors_corrected, int nb_cores, int * nbliving)
        : _bloom(bloom), _outbank(outbank), _bloocoo(bloocoo),
          _total_nb_errors_corrected (nb_errors_corrected), _local_nb_errors_corrected(0),
          _synchro(System::thread().newSynchronizer()),_temp_nb_seq_done(0), _nb_living(nbliving)
    {
        _bankwriter = new OrderedBankWriter(outbank,nb_cores*10000);
        _thread_id = __sync_fetch_and_add (_nb_living, 1);

        _tab_multivote = (unsigned char *) malloc(TAB_MULTIVOTE_SIZE*sizeof(unsigned char)); // pair of muta  = 16 nt *128 pos * 16 (max dist)
      //  printf("------- CorrectReads Custom Constructor  %p --------- tid %i \n",this,_thread_id);
        _nb_less_restrictive_correction = 3;

    }

    ~CorrectReads ()
    {
      //  printf("------- CorrectReads Destructor  %p --------- tid %i \n",this,_thread_id);

        /** We increase the global number of corrected errors. */
      //  printf("local nb errors %lli \n",_local_nb_errors_corrected);
        _thread_id = __sync_fetch_and_add (_total_nb_errors_corrected, _local_nb_errors_corrected);
        
        
#ifndef SERIAL
        _bankwriter->incDone(_temp_nb_seq_done);
        _bankwriter->waitForWriter();
#endif

        
        int nb_remaining = __sync_fetch_and_add (_nb_living, -1);

#ifndef SERIAL

        if(nb_remaining==1)
        {
             _bankwriter->FlushWriter();
        }
#endif
 
     //   free(_tab_multivote); // pourquoi plante ? lobjet correct read est il copié qq part ?
        
    }

    Bloom<kmer_type> * _bloom; // the bloom containing the solid kmers
    Bag<Sequence> *             _outbank; // the bank cto insert the result : corrected reads
    Bloocoo *         _bloocoo; // the parent bloocoo object
    unsigned char *   _tab_multivote;
    int *  _nb_living;
    int _thread_id;
    
    OrderedBankWriter * _bankwriter;
    std::vector<int> _corrected_pos;
	int _nb_less_restrictive_correction;
	
   // KmerModel           model;
   // KmerModel::Iterator itKmer;

    ISynchronizer* _synchro;
    ISynchronizer* getSynchro ()  { return _synchro; }

    u_int64_t *  _total_nb_errors_corrected;
    u_int64_t   _local_nb_errors_corrected;
    size_t _temp_nb_seq_done;
    
    
    //copy construct
    CorrectReads(const CorrectReads& cr) //called by dispatcher iterate to create N functors
    {

        //functors share smee bloom, bankwriter, bloocoo and synchronizer
        _bloom = cr._bloom;
        _outbank = cr._outbank;
        _bloocoo = cr._bloocoo;
        _synchro = cr._synchro;
        _bankwriter = cr._bankwriter;
        _nb_living = cr._nb_living;
        
        _tab_multivote = (unsigned char *) malloc(TAB_MULTIVOTE_SIZE*sizeof(unsigned char)); // pair of muta  = 16 nt *128 pos * 16 (max dist)

        _total_nb_errors_corrected = cr._total_nb_errors_corrected;
        _local_nb_errors_corrected =0;
        _temp_nb_seq_done = 0;
        _thread_id = __sync_fetch_and_add (_nb_living, 1);
       // printf("------- CorrectReads copy construct  from  %p  into %p --------- tid %i \n",&cr,this,_thread_id);

		_nb_less_restrictive_correction = cr._nb_less_restrictive_correction;
		
    }
    //assign
//    CorrectReads& operator=(const CorrectReads& cr)
//    {
//     //   printf("------- CorrectReads assign operator from  %p   ---------\n",&cr);
//        return *this;
//    }

    
};

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
Bloocoo::Bloocoo () : Tool("bloocoo"), _kmerSize(27), _inputBank (0)
{
   // _seq_num = 0;

	if(PRINT_STATS){
		for(int i=0; i<NB_CORRECTION_METHODS; i++){
			__correction_methods_successes[i] = 0;
            __correction_methods_calls[i] = 0;

		}
	}
#ifdef ERR_TAB
    pthread_mutex_init(&errtab_mutex, NULL);
#endif
    /** We add options specific to this tool. */
    getParser()->add (new OptionOneParam (STR_KMER_SIZE,                    "size of a kmer",   true));
    getParser()->add (new OptionOneParam (STR_URI_DB,                       "database",         true));   // not useful ?
    getParser()->add (new OptionOneParam (STR_KMER_SOLID,                   "solid kmers file", false));
    getParser()->add (new OptionOneParam (Bloocoo::STR_NB_MIN_VALID,        "min number of kmers to valid a correction", false,"2"));
    getParser()->add (new OptionOneParam (Bloocoo::STR_NB_VALIDATED_KMERS,  "number of kmers checked when correcting an error", false,"2"));
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void Bloocoo::execute ()
{
    /*************************************************/
    // We set some attributes (shortcuts).
    /*************************************************/
    _kmerSize           = getInput()->getInt (STR_KMER_SIZE);
    _solidFile          = getInput()->getStr (STR_KMER_SOLID);
    _nb_kmers_checked   = getInput()->getInt (STR_NB_VALIDATED_KMERS);
    _nb_min_valid = getInput()->getInt (STR_NB_MIN_VALID);
    _max_multimutation_distance = 6;
    _only_decrease_nb_min_valid = false; // false = descend les 2, recall plus faible
    
    /*************************************************/
    /** We create a bloom with inserted solid kmers. */
    /*************************************************/
    _bloom = createBloom ();
    LOCAL (_bloom);

    //iterate over initial file
    Bank inbank (getInput()->getStr(STR_URI_DB));

    /*************************************************/
    // We create a sequence iterator for the bank
    /*************************************************/
    Iterator<Sequence>* itSeq = createIterator<Sequence> (
        inbank.iterator(),
        inbank.estimateNbSequences(),
        "Iterating and correcting sequences"
    );
    LOCAL (itSeq);

    /*************************************************/
    // We create the corrected file
    /*************************************************/

    /** We get the basename from the provided URI (ie remove directory path and suffix). */
    string prefix = System::file().getBaseName (getInput()->getStr(STR_URI_DB));

    /** We set the filename as the base name + a specific suffix. */
    string fileName = prefix + string("_corrected.fasta");

    Bank outbank (fileName);

    BagCache<Sequence> *  outbank_cache = new BagCache<Sequence> (&outbank,10000);
    
    u_int64_t total_nb_errors_corrected = 0;
    int nb_corrector_threads_living;
    
#ifdef ERR_TAB
    //file with list of errors, for testing purposes
    string ferrfile = prefix + string ("_bloocoo_corr_errs.tab");

    _errfile = System::file().newFile (ferrfile, "wb");

    //file with list of errors, for testing purposes, with name of method responsible of correction
    string ferrfile2 = prefix + string ("_bloocoo_corr_errs_full.tab");
    
    _errfile_full = System::file().newFile (ferrfile2, "wb");
#endif
    
#if PRINT_LOG_MULTI
    string ferrfile3 = prefix + string ("_debug.tab");
    
    _debug = System::file().newFile (ferrfile3, "wb");
#endif
    /*************************************************/
    // We iterate over sequences and correct them
    /*************************************************/
    {
        TIME_INFO (getTimeInfo(), "sequences correction");

     //   CorrectReads fct2(*bloom, outbank, *this,total_nb_errors_corrected) ;
     //   itSeq->iterate (fct2); //
     //   
        
        nb_corrector_threads_living = 0 ;
        
#ifdef SERIAL
        setDispatcher (new SerialDispatcher());
#else
        setDispatcher (  new Dispatcher (getInput()->getInt(STR_NB_CORES)) );
#endif
        
        
        //replace outbank_cache with &outbank to remove usage of bagcache
        getDispatcher()->iterate (itSeq,  CorrectReads (_bloom, &outbank , this, &total_nb_errors_corrected,getInput()->getInt(STR_NB_CORES),&nb_corrector_threads_living),10000); // , 10000
        
       // printf("---after dispatcher iterate ---\n");
    }

    /*************************************************/
    // We gather some statistics.
    /*************************************************/
    getInfo()->add (1, "result");
    getInfo()->add (2, "nb errors corrected", "%ld", total_nb_errors_corrected);
    getInfo()->add (2, "corrected file",      fileName);
    if(PRINT_STATS){
    	printf("Correction methods stats (TwoSided, Agressive, Vote, MultiMutateVote, SideVote):\n");
		for(int i=0; i<NB_CORRECTION_METHODS; i++){
			printf("%i    /  %i \n", __correction_methods_successes[i],__correction_methods_calls[i] );
		}
	}
    
    
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
Bloom<kmer_type>* Bloocoo::createBloom ()
{
    TIME_INFO (getTimeInfo(), "fill bloom filter");

    double lg2 = log(2);
    float NBITS_PER_KMER = log (16*_kmerSize*(lg2*lg2))/(lg2*lg2);
    NBITS_PER_KMER = 12;
    printf("NBITS per kmer %f \n",NBITS_PER_KMER);
    u_int64_t solidFileSize = (System::file().getSize(_solidFile) / sizeof (kmer_type));

    u_int64_t estimatedBloomSize = (u_int64_t) ((double)solidFileSize * NBITS_PER_KMER);
    if (estimatedBloomSize ==0 ) { estimatedBloomSize = 1000; }

    /** We create the kmers iterator from the solid file. */
    Iterator<Kmer<kmer_type> >* itKmers = createIterator<Kmer<kmer_type> > (
        new IteratorFile<Kmer<kmer_type> > (_solidFile),
        solidFileSize,
        "fill bloom filter"
    );
    LOCAL (itKmers);

    /** We instantiate the bloom object. */    
    BloomBuilder<kmer_type> builder (estimatedBloomSize, 7,tools::collections::impl::BloomFactory::CacheCoherent,getInput()->getInt(STR_NB_CORES));
    Bloom<kmer_type>* bloom = builder.build (itKmers);

    /** We return the created bloom filter. */
    return bloom;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void Bloocoo::update_nb_errors_corrected(int nb_errors_corrected, u_int64_t* _local_nb_errors_corrected, bool* continue_correction){
    *_local_nb_errors_corrected += nb_errors_corrected;
    if(nb_errors_corrected > 0){
    	*continue_correction = true;
    }
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
//// twoSidedCorrection
// Method 1 for read correction
// Correct an error in an untrusted zone of size=kmerSize
//pos is the position of the error
int Bloocoo::twoSidedCorrection(int pos, char *readseq, kmer_type* kmers[], Sequence& cur_seq,std::vector<int>& tab_corrected_pos)
{
	if(PRINT_STATS){ __correction_methods_calls[TWO_SIDED] ++; }
	if(!is_pos_correctable(pos, readseq,tab_corrected_pos)){
		if(PRINT_DEBUG){ __badReadStack += "\t\t\tfailed (position is already corrected)\n"; }
		return 0;
	}
	
    kmer_type left_most_kmer = *kmers[pos-_kmerSize+1];
    kmer_type right_most_kmer = *kmers[pos];
    
    int original_nt = (readseq[pos]>>1)&3;
    int good_nt;
    kmer_type mutated_kmer, mutated_kmer_min;
    int nb_alternative = 0;

	/*
	if(pos > 65 || pos < 33){
		printf("%i\n", pos);
	}*/
	
    for(int nt=0; nt<4; nt++){
		if(nt == original_nt){
			continue;
		}
		
		//Leftmost kmer, we mutate the last NT of this kmer
		mutated_kmer = left_most_kmer;
		mutate_kmer(&mutated_kmer, 0, nt);
		mutated_kmer_min = min(mutated_kmer, revcomp(mutated_kmer, _kmerSize));
		

        
		//Check if the mutated left most kmer is trusted or not
        if(_bloom->contains(mutated_kmer_min)){
        
            //The mutated left most kmer is indexed, we can now check if the mutated right most kmer is trusted
            //Rightmost kmer, we mutate the first NT of this kmer
            mutated_kmer = right_most_kmer;
            mutate_kmer(&mutated_kmer, _kmerSize-1, nt);
            mutated_kmer_min = min(mutated_kmer, revcomp(mutated_kmer, _kmerSize));
            
            if(_bloom->contains(mutated_kmer_min)){
            	//Both the left most and right most mutated kmers are trusted so the current NT is a correction alternative
                nb_alternative += 1;
                good_nt = nt;
            }
        }
    }
    
    if(nb_alternative == 1){
    	int nb_correction = apply_correction(readseq, pos, good_nt,TWO_SIDED,cur_seq,tab_corrected_pos);
    	if(PRINT_STATS){ __correction_methods_successes[TWO_SIDED] += nb_correction; }
        return nb_correction;
    }
    
    if(PRINT_DEBUG){
    	if(nb_alternative==0){
    		__badReadStack += "\t\t\tfailed (no alternative)\n";
    	}
    	else{
    		__badReadStack += "\t\t\tfailed (multiple alternative)\n";
    	}
    }
    
    return 0;
    
}



/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
//// aggressiveCorrection
// Method 2 for read correction
// This method can correct an error in an untrusted zone of size > kmerSize (close errors)
//
// nb_kmer_check the number of additional kmers checked
// pos = position of the solid kmer next to the untrusted zone
// kmer_begin  = solid kmer at pos
// Direction:
//     RIGHT: We start from the last trusted kmer BEFORE the untrusted zone and try to correct the first error in this zone.
//     LEFT:  We start from the first trusted kmer AFTER the untrusted zone and try to correct the last error in this zone.
int Bloocoo::aggressiveCorrection(int pos, char *readseq, kmer_type* kmers[], int nb_kmer_check, int readlen, Direction direction, Sequence& cur_seq,std::vector<int>& tab_corrected_pos, int nb_kmer_offset)
{
    if(PRINT_STATS){ __correction_methods_calls[AGRESSIVE] ++; }

	//Determine corrected position depending of a left or right agressive correction
    int corrected_pos;
    if(direction ==RIGHT){
        corrected_pos = pos+_kmerSize-1;
    }
    else{
        corrected_pos = pos;
    }
    
    if(PRINT_DEBUG){ 
        std::ostringstream oss;
		oss << corrected_pos;
    	__badReadStack += "\t\t\tTry to correct pos:" + oss.str() + "\n";
    }
    

	
    //printf("%i =>    ", corrected_pos);
    
    //Cancelled the vote if the corrected position is not correctable
	if(!is_pos_correctable(corrected_pos, readseq,tab_corrected_pos)){
		if(PRINT_DEBUG){ __badReadStack += "\t\t\tfailed (position is already corrected)\n"; }
		return 0;
	}
	
	/*
    if(corrected_pos==0){
    	if(PRINT_DEBUG){__badReadStack += "\t\tErr pos 0 => Read Side correction\n"; }
    	//printf("start read side correction (err_pos: %i)\n", corrected_pos);
    	//printf("seq_num: %i\n", _seq_num);
    	int nb_cor = extendedAgressiveCorrection(pos, readseq, kmers, 2, LEFT,cur_seq,tab_corrected_pos);
    	return nb_cor;
    }
    else if(corrected_pos==readlen-1){
    	if(PRINT_DEBUG){__badReadStack += "\t\tErr pos 99 => Read Side correction\n"; }
    	//printf("start read side correction (err_pos: %i)\n", corrected_pos);
    	//printf("seq_num: %i\n", _seq_num);
    	int nb_cor = extendedAgressiveCorrection(pos, readseq, kmers, 2, RIGHT,cur_seq,tab_corrected_pos);
    	return nb_cor;
    }*/
    
    //Determine the minimum vote threshold
    //It's always the param _nb_min_valid except for the read sides (Sides have limited number of kmers to check)
    //int min_valid = _nb_min_valid-nb_kmer_offset;
    
    int current_max_nb_checkable;
    if(direction ==RIGHT){
    	current_max_nb_checkable = readlen - corrected_pos;
    }
    else{
    	current_max_nb_checkable = corrected_pos + 1;
    }
    
    int vote_threshold = _nb_min_valid-nb_kmer_offset;
    if(!_only_decrease_nb_min_valid) nb_kmer_check -= nb_kmer_offset;
    nb_kmer_check = max(1, nb_kmer_check);
    vote_threshold = max(1, vote_threshold);
    
    //printf("Pos: %i     nb_kmer_check: %i     max_kmer_checkable: %i     vote_treshold: %i\n", corrected_pos, nb_kmer_check, current_max_nb_checkable, vote_threshold);
    
    //int nb_extra_kmer = max(0, nb_kmer_check-current_max_nb_checkable);
    //nb_extra_kmer = min(1, nb_extra_kmer);
    int nb_extra_kmer = max(0, vote_threshold-current_max_nb_checkable);
    //nb_extra_kmer = min(2, nb_extra_kmer);
    
    //current_max_nb_checkable += nb_extra_kmer;
    //current_max_nb_checkable += nb_extra_kmer;
    //current_max_nb_checkable = nb_extra_kmer + nb_kmer_check;
    
    
    nb_kmer_check = min(current_max_nb_checkable, nb_kmer_check);
    
    //int vote_threshold = min(nb_extra_kmer + nb_kmer_check, _nb_min_valid-nb_kmer_offset);
    //nb_kmer_check = max(1, nb_kmer_check);
    
    //printf("Pos: %i      nb_extra_kmer: %i    nb_kmer_check: %i        vote_treshold: %i\n\n", corrected_pos, nb_extra_kmer, nb_kmer_check, vote_threshold);
    
    //printf("%i %i %i\n", corrected_pos, vote_threshold, nb_kmer_check);
    
    
    //If the number of checkable kmers is inferior to the vote threshold then the vote is cancelled
	if(nb_kmer_check+nb_extra_kmer < vote_threshold){
		if(PRINT_DEBUG){
            std::ostringstream oss;
			oss << nb_kmer_check;
			std::ostringstream oss2;
			oss2 << vote_threshold;
            __badReadStack += "\t\t\tfailed (nb_kmer_checked < nb_min_valid) " + oss.str() + " < " + oss2.str() + " \n";
        }
		return 0;
	}
	
	
	
	//printf("%i    nb_checked %i,    threshold %i\n", corrected_pos, nb_kmer_check, vote_threshold);
	
	kmer_type kmer_begin = *kmers[pos];
	int votes [4] = {0, 0, 0, 0};
    int good_nt;
    int nb_possibles_firstk=0;
	char nt_temp;
    
    KmerModel model (_kmerSize);
    kmer_type current_kmer, current_kmer_min;
    
    int original_nt;
    if(direction ==RIGHT){
        original_nt = (readseq[corrected_pos]>>1)&3;
    }
    else{
    	original_nt = (readseq[corrected_pos]>>1)&3;
    }
        
        
    /*
    int pre_good_nt = -1;
    int pre_nb_alternative = 0;
      
    for(int nt=0; nt<4; nt++)
    {
		if(nt == original_nt){
			continue;
		}
			
		current_kmer = kmer_begin;
		
		if(direction ==RIGHT){
			mutate_kmer(&current_kmer, 0, nt);
		}
		else{
		    mutate_kmer(&current_kmer, _kmerSize-1, nt);
		}
		
		//current_kmer = codeSeedBin(&model, &kmer_begin, nt, direction);
        current_kmer_min = min(current_kmer, revcomp(current_kmer, _kmerSize));
        
        if(_bloom->contains(current_kmer_min)) //kmer is indexed
        {
			pre_nb_alternative += 1;
			pre_good_nt = nt;
			//votes[nt] += 1;
            //nb_possibles_firstk++;
        }
        else{
        	continue; // if first kmer is not valid, this possible nt is eliminated
        }
	}
	
	if(pre_nb_alternative == 1){
		int nb_correction = apply_correction(readseq, corrected_pos, pre_good_nt, AGRESSIVE,cur_seq,tab_corrected_pos);
    	if(PRINT_STATS){ __correction_methods_successes[AGRESSIVE] += nb_correction; }
        return nb_correction;
	}*/
	
	
	
    for(int nt=0; nt<4; nt++)
    {
		if(nt == original_nt){
			continue;
		}
			
		current_kmer = kmer_begin;
		
		if(direction ==RIGHT){
			mutate_kmer(&current_kmer, 0, nt);
		}
		else{
		    mutate_kmer(&current_kmer, _kmerSize-1, nt);
		}
		
		//current_kmer = codeSeedBin(&model, &kmer_begin, nt, direction);
        current_kmer_min = min(current_kmer, revcomp(current_kmer, _kmerSize));
        
        if(_bloom->contains(current_kmer_min)) //kmer is indexed
        {
			votes[nt] += 1;
            nb_possibles_firstk++;
        }
        else{
        	//continue; // if first kmer is not valid, this possible nt is eliminated
        }

		//nb_kmer_check-1: -1 because the first mutation above is counted as a kmer check
		for (int ii=0; ii<nb_kmer_check-1; ii++)
		{
		    
		    if(direction == RIGHT){
		        nt_temp = readseq[corrected_pos+1+ii];
		    }
		    else{
		        nt_temp = readseq[corrected_pos-1-ii];
		    }
		    
		    codeSeedNT(&model, &current_kmer, nt_temp, direction);
		    
		    current_kmer_min = min(current_kmer, revcomp(current_kmer, _kmerSize));
		    if(_bloom->contains(current_kmer_min)) //kmer is indexed
		    {
		        votes[nt] += 1;
		    }
		    
		}
		
		if(nb_extra_kmer > 0 /*&& votes[nt] > 0*/){
			extendedAgressiveCorrection(votes, &model, &current_kmer, nt, nb_extra_kmer, direction);
		}
    }


    //print_agressive_votes(votes);

	//Searching max score in votes
	int max_score = votes[0];
	for(int i=1; i<4; i++){
 		if(votes[i] > max_score){
			max_score = votes[i];
		}
	}

	//int t = (nb_kmer_check*25)/100;
	//int t = max(1, nb_kmer_check);
	//int t = 1; //max(1, nb_kmer_check);
	//if(nb_possibles_firstk==1) vote_threshold=min(current_max_nb_checkable,3);
    
	if(max_score < vote_threshold){
		if(PRINT_DEBUG){ __badReadStack += "\t\t\tfailed (max score %i is too low)\n"; }
		return 0;
	}
	
	//printf("max: %i\n", max_score);

	//We have the max score in votes, now we have to check if this score
	//is uniq or if more than one nt has the max score
	int nb_alternative = 0;

	for(int i=0; i<4; i++){
 		if(votes[i] == max_score){
			nb_alternative += 1;
			good_nt = i;
		}
	}

	//printf("good_nt: %i, nb_alternative: %i\n", good_nt, nb_alternative);

	//if max score is uniq we can correct the sequence
	//else if nb_nt_max > 1 then there are more than one candidate for the correction
	//so we cannot determine a good correction, the correction is cancelled.
	if(nb_alternative == 1){
		int nb_correction = apply_correction(readseq, corrected_pos, good_nt,AGRESSIVE,cur_seq,tab_corrected_pos);
    	if(PRINT_STATS){ __correction_methods_successes[AGRESSIVE] += nb_correction; }
        return nb_correction;
	}

    
    if(PRINT_DEBUG){ __badReadStack += "\t\t\tfailed (multiple alternative)\n"; }
    
    return 0;
    
}



/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
int Bloocoo::voteCorrectionInUntrustedZone(int start_pos, int end_pos, char *readseq, kmer_type* kmers[], int nb_kmer_checked, Sequence& cur_seq,std::vector<int>& tab_corrected_pos, int nb_kmer_offset){
	//printf("%i:    %i %i\n", _seq_num, start_pos, end_pos);
	//printf("\n\t");
	//(*kmers[start_pos]).printASCII(_kmerSize);
	//(*kmers[end_pos]).printASCII(_kmerSize);
	
	//start_pos = max(0, start_pos);
	
	int nb_errors_cor = 0;
    if(PRINT_STATS){ __correction_methods_calls[VOTE] ++; }

	while(nb_errors_cor == 0 && start_pos < end_pos){
		int untrusted_zone_size = end_pos - start_pos;
		int new_nb_checked = min(nb_kmer_checked, untrusted_zone_size);
		if(new_nb_checked <= 0){ break; }
		
		if(PRINT_DEBUG){
			std::ostringstream oss;
			oss << start_pos;
			std::ostringstream oss2;
			oss2 << end_pos;
		 	__badReadStack += "\tVote pos: " + oss.str() + " " + oss2.str() + "\n";
		}
		
		nb_errors_cor = voteCorrection(start_pos, readseq, kmers, new_nb_checked,cur_seq,tab_corrected_pos, nb_kmer_offset);
		start_pos += _kmerSize/2;//new_nb_checked;
	}
	
	return nb_errors_cor;
}
		            		
/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
//// voteCorrection
// Method 3 for read correction
// This method is used if the method 2 fails to correct an error
//
// nb_kmer_check the number of additional kmers checked
// pos: the last trusted position in readseq before the untrusted zone
// kmer_begin: solid kmer at pos, it's the last trusted kmer before the untrusted zone
int Bloocoo::voteCorrection(int start_pos, char *readseq, kmer_type* kmers[], int nb_kmer_check, Sequence& cur_seq,std::vector<int>& tab_corrected_pos, int nb_kmer_offset)
{
	int vote_threshold = _nb_min_valid-nb_kmer_offset;
	if(!_only_decrease_nb_min_valid) nb_kmer_check -= nb_kmer_offset;
	nb_kmer_check = max(1, nb_kmer_check);
	vote_threshold = max(1, vote_threshold);
	
	if(nb_kmer_check < vote_threshold){
		if(PRINT_DEBUG){ __badReadStack += "\t\t\tfailed (nb_kmer_checked < nb_min_valid)\n";}
		return 0;
	}
	
    int good_nt;
    //KmerModel model (_kmerSize);
    //kmer_type kmer_begin = *kmers[start_pos];
    kmer_type current_kmer, current_kmer_min, mutated_kmer;

	int nb_column = nb_kmer_check + _kmerSize - 1;
	int votes[nb_column][4];
	for (int i = 0; i < nb_column; i++) {
		for (int nt=0; nt < 4; nt++) {
			votes[i][nt] = 0;
		}
	}
	
	//kmer_begin.printASCII(_kmerSize);
	//current_kmer = kmer_begin;
	
	for(int i=0; i < nb_kmer_check; i++){

		int read_pos = start_pos + i;
		
		current_kmer = *kmers[read_pos];
		//current_kmer = model.codeSeedRight(current_kmer, readseq[read_pos+(_kmerSize-1)], Data::ASCII, KMER_DIRECT);
		
		//printf("\tcurrent kmer:    "); current_kmer.printASCII(_kmerSize);
		
		current_kmer_min = min(current_kmer, revcomp(current_kmer, _kmerSize));
		//if(_bloom->contains(current_kmer_min)){
		//    continue;
		//}
		//if(!is_pos_correctable(read_pos, readseq)){
		//	continue;
		//}
				
		for(int kpos=0; kpos < _kmerSize; kpos++){
			if(!is_pos_correctable(read_pos+kpos, readseq,tab_corrected_pos)){
				continue;
			}
		
			int original_nt = (readseq[read_pos+kpos]>>1)&3;
			
    		for(int nt=0; nt<4; nt++){
				if(nt == original_nt){
					continue;
				}
				
				mutated_kmer = current_kmer;
				mutate_kmer(&mutated_kmer, _kmerSize-1-kpos, nt);
				//mutated_kmer.printASCII(_kmerSize);
				
				/*
				printf("mutate at %i with %c\n", _kmerSize-1-kpos, bin2NT[nt]);
				printf("curre_kmer: ");
				current_kmer.printASCII(_kmerSize);
				printf("mutated_kmer: ");
				mutated_kmer.printASCII(_kmerSize);
				printf("\n");
				*/
				
				current_kmer_min = min(mutated_kmer, revcomp(mutated_kmer, _kmerSize));
				if(_bloom->contains(current_kmer_min))
				{
				    votes[i+kpos][nt] += 1;
					//print_votes(votes);
				}
    		}
		}
	}

	//print_votes(votes, nb_column);

	//Search max vote in the matrix
	int vote, maxVote = 0 , nb_max_vote= 0;
	for (int i = 0; i < nb_column; i++) {
		for (int nt=0; nt < 4; nt++) {
			vote = votes[i][nt];
            if (vote == maxVote) {
                nb_max_vote ++ ;
			}
			if (vote > maxVote) {
				maxVote = vote;
                nb_max_vote =1;
			}
		}
	}

    
    if(nb_max_vote !=1 ){
		if(PRINT_DEBUG){ __badReadStack += "\t\t\tfailed nb_max_vote !=1 \n"; }
		return 0;
	}
    
	//int t = (nb_kmer_check*25)/100;
	//t = max(2, t);
	//int t = max(1, nb_kmer_check);
	//int t = _nb_min_valid;
	if(maxVote < vote_threshold){
		if(PRINT_DEBUG){ __badReadStack += "\t\t\tfailed\n"; }
		return 0;
	}
	
	//printf("max vote: %i", maxVote);
	int nb_alternative, corrected_pos;
	int nb_correction = 0;

	
    //si plusieurs max vote, on les corrige tous : ptet cest pas bien : ils sont peut etre plus proches qun kmer, normalement impossible avec  1 muta
	
	//For each column in the matrix, if the max vote is find and is uniq we can apply the correction
	for (int i = 0; i < nb_column; i++) {
		nb_alternative = 0;
		for (int nt=0; nt < 4; nt++) {
			vote = votes[i][nt];
			if (vote == maxVote) {
				nb_alternative += 1;
				good_nt = nt;
			}
		}
		if(nb_alternative == 1){
			corrected_pos = start_pos + i;
			int nb_cor = apply_correction(readseq, corrected_pos, good_nt,VOTE,cur_seq,tab_corrected_pos);
			nb_correction += nb_cor;
	    	if(PRINT_STATS){ __correction_methods_successes[VOTE] += nb_cor; }
		}
	}

	if(PRINT_DEBUG){
		if(nb_correction == 0){
			__badReadStack += "\t\t\tfailed (multiple alternative)\n";
		}
	}
	
    return nb_correction;
    
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
int Bloocoo::apply_correction(char *readseq, int pos, int good_nt,int algo, Sequence& cur_seq,std::vector<int>& tab_corrected_pos){
	if(!is_pos_correctable(pos, readseq,tab_corrected_pos)){
		if(PRINT_DEBUG){ __badReadStack += "\t\t\tFailed in apply_correction (non correctable pos)\n"; }
		return 0;
	}
	
    static const char* methodname[] = { "twosided", "aggressive", "vote", "multi" , "side"  };

	if(PRINT_DEBUG){
		std::ostringstream oss;
		oss << cur_seq.getIndex();
		std::string str__seq_num = oss.str();
		std::ostringstream oss2;
		oss2 << pos;
		std::string str_pos = oss2.str();
		__badReadStack += "\t\t\t" + str__seq_num + "\t" + str_pos + "\t" + bin2NT[good_nt] + "\t" + readseq[pos] + "\n";
	}
	//printf("\t\t%i\t%i\t%c\t%c\n",_seq_num,pos,bin2NT[good_nt],readseq[pos]);
	
#ifdef ERR_TAB
    pthread_mutex_lock(&errtab_mutex);

    _errfile->print("%i\t%i\t%c\t%c\n",cur_seq.getIndex(),pos,bin2NT[good_nt],readseq[pos]);
    _errfile_full->print("%i\t%i\t%c\t%c;%s\n",cur_seq.getIndex(),pos,bin2NT[good_nt],readseq[pos],methodname[algo]);
    pthread_mutex_unlock(&errtab_mutex);
#endif
    
    readseq[pos] = bin2NT[good_nt]; //correc sequence in ram
    // printf("error found pos %i  read %i\n",ii, _bloocoo->_seq_num);
    tab_corrected_pos.push_back(pos);
    return 1;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void Bloocoo::codeSeedBin(KmerModel* model, kmer_type* kmer, int nt, Direction direction){
	if(direction == RIGHT){
		*kmer = model->codeSeedRight(*kmer, nt, Data::INTEGER, KMER_DIRECT);
	}
	else{
		kmer_type kmer_end_rev = revcomp(*kmer, _kmerSize);
		*kmer = model->codeSeedRight (kmer_end_rev, binrev[nt], Data::INTEGER, KMER_DIRECT);
		*kmer = revcomp(*kmer, _kmerSize);
	}
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void Bloocoo::codeSeedNT(KmerModel* model, kmer_type* kmer, char nt, Direction direction){
	//Attention ici, penser à utiliser une fonction binToNT statique.
	return codeSeedBin(model, kmer, (nt>>1)&3, direction);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
//// print_agressive_votes
//Debug function for agressive correction, print the votes tab
void Bloocoo::print_agressive_votes(int votes[4]){
	printf("#############################################\n");
	for(int i=0; i<4; i++){
		printf("%c: %i\n", bin2NT[i], votes[i]);
	}
	printf("#############################################\n");
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
//print_votes
//Debug function for vote correction, print the votes matrix
void Bloocoo::print_votes(int votes[][4], int nb_column){
	printf("#############################################\n");
	for(int nt=0; nt<4; nt++){
		printf("%c  ", bin2NT[nt]);
		for(int i=0; i<nb_column; i++){
			printf("%i ", votes[i][nt]);
		}
		printf("\n");
	}
	printf("#############################################\n");
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void Bloocoo::print_read_correction_state(KmerModel* model, Sequence& s){
	KmerModel::Iterator itKmer2 (*model);
	kmer_type kmer, kmer_min;
	
	itKmer2.setData (s.getData(),KMER_DIRECT);
	int i3=0;
	//printf("\nread correction state (%i):\n", seq_num);
	//printf("\t%s\n", s.getDataBuffer());
	//printf("\t");
	
	std::ostringstream oss;
    oss << s.getIndex();


	__badReadStack += "Read correction state (" + oss.str() += "):\n";
	__badReadStack += "\t" + std::string(s.getDataBuffer());
	__badReadStack += "\n\t";
	
	for (itKmer2.first(); !itKmer2.isDone(); itKmer2.next(),i3++){
		kmer = *itKmer2;
		kmer_min = std::min(  revcomp (kmer, _kmerSize),kmer );
		if (_bloom->contains(kmer_min)){
			//printf("1");
			__badReadStack += "1";
		}
		else{
			//printf("0");
			__badReadStack += "0";
		}
	}
	//printf("\n");
	__badReadStack += "\n";
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void Bloocoo::print_read_if_not_fully_corrected(KmerModel* model, Sequence& s){
	std::string read_state = "";
	std::string normal_state = "";
	KmerModel::Iterator itKmer2 (*model);
	kmer_type kmer, kmer_min;
	
	itKmer2.setData (s.getData(),KMER_DIRECT);
	int i3=0;
	//printf("\nread correction state (%i):\n", _seq_num);
	//printf("\t%s\n", s.getDataBuffer());
	//printf("\t");
	for (itKmer2.first(); !itKmer2.isDone(); itKmer2.next(),i3++){
		normal_state += "1";
		kmer = *itKmer2;
		kmer_min = std::min(  revcomp (kmer, _kmerSize),kmer );
		if (_bloom->contains(kmer_min)){
			read_state += "1";
		}
		else{
			read_state += "0";
		}
	}
	if(normal_state.compare(read_state) != 0){
		//print_read_correction_state(_bloom, model, _kmerSize, s, _seq_num);
		printf("%s", __badReadStack.c_str());
	}
	
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
//// mutate_kmer
// fonction to mutate kmer : takes kmer, pos (de la fin ,0-based), and nt (nt = 0,1,2ou 3) 
// par exemple  kmer ,2 , C    avec kmer =  AAAAAAAAAA
//
// return : 
// AAAAAAACAA
void Bloocoo::mutate_kmer(kmer_type * kmer, int pos, char nt)
{
    kmer_type reset_mask =  ~((kmer_type)3 << (pos*2));
    kmer_type set_mask =  ((kmer_type)nt) << (pos*2);
    *kmer = (*kmer & reset_mask )  | set_mask;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
bool Bloocoo::is_pos_correctable(int pos, char* readseq,std::vector<int>& tab_corrected_pos){
	bool N_at_pos = (readseq[pos] == 'N');
	bool is_pos_corrected = (std::find(tab_corrected_pos.begin(), tab_corrected_pos.end(), pos) != tab_corrected_pos.end());
	return !N_at_pos && !is_pos_corrected;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
int Bloocoo::multiMutateVoteCorrectionInUntrustedZone(int start_pos, int end_pos, char *readseq, kmer_type* kmers[], int nb_kmer_checked, unsigned char* _tab_multivote, int expected_first_pos, int expexted_second_pos, int readlen, Sequence& cur_seq,std::vector<int>& tab_corrected_pos, int nb_kmer_offset){
	//printf("%i:    %i %i\n", _seq_num, start_pos, end_pos);
	//printf("\n\t");
	//(*kmers[start_pos]).printASCII(_kmerSize);
	//(*kmers[end_pos]).printASCII(_kmerSize);
	
	//start_pos = max(0, start_pos);
	//printf("\n\tMulti mutate start !!!!!!!!!!\n");
	
	//_tab_multivote = (unsigned char *) malloc(TAB_MULTIVOTE_SIZE*sizeof(unsigned char)); // pair of muta  = 16 nt *128 pos * 16 (max dist)
    if(PRINT_STATS){ __correction_methods_calls[MULTI_MUTATE_VOTE] ++; }

#if PRINT_LOG_MULTI
    _debug->print("%i\t",end_pos-start_pos +1);
#endif
	int nb_errors_cor = 0;
	int nz=0;
	while(nb_errors_cor == 0 && start_pos < end_pos){
		int untrusted_zone_size = end_pos - start_pos;
		int new_nb_checked = min(nb_kmer_checked, untrusted_zone_size);
		if(new_nb_checked <= 0){ break; }
		
		if(PRINT_DEBUG){
			std::ostringstream oss;
			oss	<< start_pos;
			std::ostringstream oss2;
			oss2 << end_pos;
		 	__badReadStack += "\tMulti Mutate Vote pos: " + oss.str() + " " + oss2.str() + "\n";
		}
		
        nz++;
		nb_errors_cor = multiMutateVoteCorrection(start_pos, readseq, kmers, new_nb_checked, _tab_multivote,expected_first_pos,expexted_second_pos,readlen,cur_seq,tab_corrected_pos, nb_kmer_offset);
		start_pos += new_nb_checked;
		// start_pos += _kmerSize/2;
	}
#if PRINT_LOG_MULTI
    _debug->print("%i\t%i\n",nz,nb_errors_cor);
#endif
	return nb_errors_cor;
}


/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
int Bloocoo::multiMutateVoteCorrection(int start_pos, char *readseq, kmer_type* kmers[], int nb_kmer_check, unsigned char* _tab_multivote,int expected_first_pos, int expected_second_pos, int readlen, Sequence& cur_seq,std::vector<int>& tab_corrected_pos, int nb_kmer_offset)
{
    
    int current_max_nb_checkable = 999999;
    if(expected_second_pos > 0 &&  expected_second_pos < 10){
        current_max_nb_checkable = 3;
    }
    if(expected_first_pos > 0   &&  (readlen -expected_first_pos)<10 ){
    	current_max_nb_checkable = 3;
    }
    
    int vote_threshold = min(current_max_nb_checkable, _nb_min_valid-nb_kmer_offset);

	//int vote_threshold = _nb_min_valid-kmer_offset;
	if(!_only_decrease_nb_min_valid) nb_kmer_check -= nb_kmer_offset;
	nb_kmer_check = max(1, nb_kmer_check);
	vote_threshold = max(1, vote_threshold);
	

	if(nb_kmer_check < vote_threshold){
		if(PRINT_DEBUG){ __badReadStack += "\t\t\tfailed (nb_kmer_checked < nb_min_valid)\n";}
		return 0;
	}
	
	memset(_tab_multivote, 0, TAB_MULTIVOTE_SIZE*sizeof(unsigned char)); 
    
	int nb_max_vote = 0;
	int max_vote = 0;
	int vote;
	int nb_correction = 0;
	int good_index;
	
	//int nb_column = nb_kmer_check + _kmerSize - 1;
	//std::map<std::string, int> votes;
	
	for(int i=0; i < nb_kmer_check; i++){
		kmer_type current_kmer = *kmers[start_pos+i];
		multiMutateVoteCorrectionRec(start_pos, i, current_kmer, readseq, kmers, nb_kmer_check, 0, 0, _tab_multivote, &max_vote, &nb_max_vote, &good_index, 0, 0,expected_first_pos,expected_second_pos,cur_seq,tab_corrected_pos);
		//multiMutateVoteCorrectionRec(start_pos, i, current_kmer, readseq, kmers, nb_kmer_check, 0, 0, _tab_multivote, 0, 0);
	}

    
    
	//Search max vote in the hash

	/*
	for(int i=0; i<TAB_MULTIVOTE_SIZE; i++){
		vote = _tab_multivote[i];
		if (vote > max_vote) {
			max_vote = vote;
		}
		//if(_tab_multivote[i] >= 1){
			
			//decode_index(_tab_multivote[i], int * pos1, int * dist, int * nt1, int * nt2)
			//printf("%i", _tab_multivote[i]);
		//}
	}*/
	
	if(max_vote < vote_threshold){
		if(PRINT_DEBUG){ __badReadStack += "\t\t\tfailed (max vote too low)\n"; }
		return 0;
	}
	
	/*
	for(int i=0; i<TAB_MULTIVOTE_SIZE; i++){
		vote = _tab_multivote[i];
		if (vote == max_vote) {
			nb_max_vote += 1;
			good_index = i;
		}
	}*/
	
	if(nb_max_vote != 1){
		if(PRINT_DEBUG){ __badReadStack += "\t\t\tfailed max vote is not uniq\n"; }
		return 0;
	}
	

	
	

	
	
	
	//printf("%i\n", nb_max_vote);

	
	
	//printf("--------------- %i \n", maxVote);
	int pos1, dist, nt1, nt2;
	decode_index(good_index, &pos1, &dist, &nt1, &nt2);
	//printf("\t%i %c %i %c\n", pos1, bin2NT[nt1], dist, bin2NT[nt2]);
	
	int nb_cor = apply_correction(readseq, start_pos+pos1, nt1,MULTI_MUTATE_VOTE,cur_seq,tab_corrected_pos);
	if(PRINT_STATS){ __correction_methods_successes[MULTI_MUTATE_VOTE] += nb_cor; }
	nb_correction += nb_cor;
		
	nb_cor = apply_correction(readseq, start_pos+pos1+dist, nt2,MULTI_MUTATE_VOTE,cur_seq,tab_corrected_pos);
	if(PRINT_STATS){ __correction_methods_successes[MULTI_MUTATE_VOTE] += nb_cor; }
	nb_correction += nb_cor;
		
	//printf("\t%i %c %i %c\n", start_pos+pos1, bin2NT[nt1], start_pos+pos1+dist, bin2NT[nt2]);
	
	/*
	std::string good_mutations;
	
	std::map<std::string, int>::iterator iter;
	
	for (iter = votes.begin(); iter != votes.end(); ++iter) {
		//printf("%s\n", iter->first.c_str());
		int vote = iter->second;
		if (vote > maxVote) {
			maxVote = vote;
		}
	}
	
	//int t = max(3, nb_kmer_check-1);
	//int t = _nb_min_valid;
	if(maxVote < vote_threshold){
		if(PRINT_DEBUG){ __badReadStack += "\t\t\tfailed\n"; }
		return 0;
	}
	
	//check if max is uniq !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!(Optimisation possible)	
	for (iter = votes.begin(); iter != votes.end(); ++iter) {
		//printf("%s\n", iter->first.c_str());
		int vote = iter->second;
		if (vote == maxVote) {
			nbMax += 1;
			good_mutations = iter->first;
		}
	}
		
	if(nbMax != 1){
		return 0;
	}
	
	//Parsing the good mutation chain
	//Mutations chain has the following pattern ("pos1 nt1 pos2 nt2 pos3 nt3...")
	std::string token;
	//printf("---------- %s\n", good_mutations.c_str());
	std::stringstream ss(good_mutations);
	int pos, nt;
	int nb_correction = 0;
	int i=0;
	
	while (ss >> token){
		//printf("%s\n", token.c_str());
		i += 1;
		
		if(i==1){ //Storing mutation position
			_iss.clear();
			_iss.str(token);
			_iss >> pos;
		}
		else if(i==2){ //Storing mutation nt and applied correction
			_iss.clear();
			_iss.str(token);
			_iss >> nt;
			int nb_cor = apply_correction(readseq, start_pos+pos, nt);
			if(PRINT_STATS){ __correction_methods_successes[MULTI_MUTATE_VOTE] += nb_cor; }
			nb_correction += nb_cor;
			i = 0;
		}
		
	}
	*/
	return nb_correction;
    
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
//int Bloocoo::multiMutateVoteCorrectionRec(int start_pos, int kmer_offset, kmer_type current_kmer, char *readseq, kmer_type* kmers[], int nb_kmer_check, int kmer_index, int current_nb_mutation, unsigned char* _tab_multivote, int* max_vote, int* nb_max_vote, int *good_index, int pos1, int nt1){
void Bloocoo::multiMutateVoteCorrectionRec(int start_pos, int kmer_offset, kmer_type current_kmer, char *readseq, kmer_type* kmers[], int nb_kmer_check, int kmer_index, int current_nb_mutation, unsigned char* _tab_multivote, int* max_vote, int* nb_max_vote, int *good_index, int pos1, int nt1, int expected_first_pos, int expected_second_pos, Sequence& cur_seq,std::vector<int>& tab_corrected_pos){
	int read_pos = start_pos + kmer_offset;
	
	
	//Determine the distance between 2 mutations
	int end_kpos = _kmerSize;
	if(current_nb_mutation == 1){
		end_kpos = std::min(kmer_index+_max_multimutation_distance, (int)_kmerSize);
		//printf("%i %i\n", kmer_index, end_kpos);
	}
	
	//Iterating on kmer positions from kmer_index to the max mutate distance
	//If the max number of mutation is not reached, we iterate over all the kmer positions
	for(int kpos=kmer_index; kpos < end_kpos; kpos++){
		
		//Don't vote if the current position is not correctable
		if(!is_pos_correctable(read_pos+kpos, readseq,tab_corrected_pos)){
			continue;
		}
        
        //si celui ci a abs 1
        //expected first pos is known, restrict search to +-1
		if((expected_first_pos >=0) &&
           (current_nb_mutation == 0) &&
           abs(read_pos+kpos - expected_first_pos) > 0 ){ //only test pos +-1 of expected error pos
			continue;
		}
        
        //ceux là (extension vers left) ne semblent pas impacter recall, et gain perf correct
        //expected second pos is known  restrict first to max_dist away
        if((expected_second_pos >=0) &&
           (current_nb_mutation == 0) &&
           abs(read_pos+kpos - expected_second_pos) >  (1+_max_multimutation_distance) ){
			continue;
		}
    
        if((expected_second_pos >=0) &&
           (current_nb_mutation == 1) &&
           abs(read_pos+kpos - expected_second_pos) > 0 ){
			continue;
		}
        
        
		int original_nt = (readseq[read_pos+kpos]>>1)&3;
		
		//Iterating A, C, G, T except the original_nt
		for(int nt=0; nt<4; nt++){
			if(nt == original_nt){
				continue;
			}
			
			//Mutate the current kmer 
			kmer_type mutated_kmer = current_kmer;
			mutate_kmer(&mutated_kmer, _kmerSize-1-kpos, nt);
			
			/*
			//Add this new mutation ("pos nt") to the string that stores the mutations chain. (Ex: "60 A 12 B...")
			//_oss.clear();
			_oss.str("");
			_oss << kmer_offset+kpos;
			std::string mutations_copy = mutations + " " + _oss.str();
			//_oss.clear();
			_oss.str("");
			_oss << nt;
			mutations_copy += " " + _oss.str();
			*/
			
			//If the max number of mutations is reached, we can start voting for the mutated kmers
			if(current_nb_mutation == 1){
				/*
				printf("kmer index: %i\n", kpos);
				printf("\tcurent kmer:  "); current_kmer.printASCII(_kmerSize);
				printf("\tmutated kmer:  "); mutated_kmer.printASCII(_kmerSize);
				printf("\n");*/
				
				kmer_type mutated_kmer_min = min(mutated_kmer, revcomp(mutated_kmer, _kmerSize));
				if(_bloom->contains(mutated_kmer_min)){
					//int dist = kpos - kmer_index +1;
					//printf("%i %c %i %c\n", pos1+start_pos, bin2NT[nt1], dist, bin2NT[nt]);
					int index = make_index(pos1, kpos-kmer_index+1, nt1, nt);
					_tab_multivote[index] += 1;
					
					if(_tab_multivote[index] == *max_vote){
						*nb_max_vote += 1;
					}
					else if(_tab_multivote[index] > *max_vote){
						*max_vote = _tab_multivote[index];
						*nb_max_vote = 1;
						*good_index = index;
					}
					

				}
			}
			//If the max number of mutations is not reached, we have to recurcivelly continue to applied mutations before voting
			//The kmer_index and the current number of mutation is incremented 
			else{
				//multiMutateVoteCorrectionRec(start_pos, kmer_offset, mutated_kmer, readseq, kmers, nb_kmer_check, kpos+kmer_index+1, current_nb_mutation+1, _tab_multivote, kmer_offset+kpos, nt);
				multiMutateVoteCorrectionRec(start_pos, kmer_offset, mutated_kmer, readseq, kmers, nb_kmer_check, kpos+kmer_index+1, current_nb_mutation+1, _tab_multivote, max_vote, nb_max_vote, good_index, kmer_offset+kpos, nt,expected_first_pos,expected_second_pos,cur_seq,tab_corrected_pos);
			}
		}
	}
}




/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void Bloocoo::extendedAgressiveCorrection(int votes[4], KmerModel* model, kmer_type* last_mutated_kmer, int mutated_nt, int nb_kmer_check, Direction direction)
{
    if(PRINT_STATS){ __correction_methods_calls[READ_SIDE] ++; }

	bool posVoted[nb_kmer_check];
	
	for(int i=0; i<nb_kmer_check; i++){
		posVoted[i] = false;
	}
	
	//printf("%i\n", posHasVote[0]);
	//for(int i=0; i<nb_kmer_check; i++){
		
	//printf("------------------------------------------------\n");
	//printf("last kmer: "); (*last_mutated_kmer).printASCII(_kmerSize);
	//print_agressive_votes(votes);
	extendedAgressiveCorrectionRec(votes, model, last_mutated_kmer, mutated_nt, nb_kmer_check, direction, posVoted, 0);
	//print_agressive_votes(votes);
	//printf("------------------------------------------------\n");
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void Bloocoo::extendedAgressiveCorrectionRec(int votes[4], KmerModel* model, kmer_type* mutated_kmer, int mutated_nt, int nb_kmer_check, Direction direction, bool posVoted[], int depth)
{
	if(depth >= nb_kmer_check){
		return;
	}
	
	for(int nt=0; nt<4; nt++){
		
		kmer_type current_kmer = *mutated_kmer;
		codeSeedBin(model, &current_kmer, nt, direction);
		
		
		
		kmer_type current_kmer_min = min(current_kmer, revcomp(current_kmer, _kmerSize));
		if(_bloom->contains(current_kmer_min)) //kmer is indexed
		{
			if(!posVoted[depth]){
				//current_kmer.printASCII(_kmerSize);
				votes[mutated_nt] += 1;
				posVoted[depth] = true;
			}
			
			extendedAgressiveCorrectionRec(votes, model, &current_kmer, mutated_nt, nb_kmer_check, direction, posVoted, depth+1);
			
		}
	}
		
}



/*********************************************************************
 ** METHOD  :make_indice
 ** PURPOSE :encode mutation pair info (pos1, dist_to_pos2, nt1,nt2) into a single int
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
unsigned int Bloocoo::make_index(int pos1, int dist, int nt1, int nt2)
{
    unsigned int idx = pos1;
    
    idx <<= 4;
    idx |=  (dist & 15); //dist is constrained to [0-15]
    idx <<= 2;
    idx |=  nt1; 
    idx <<= 2;
    idx |=  nt2;
    
    return idx;
}


/*********************************************************************
 ** METHOD  :make_indice
 ** PURPOSE :encode mutation pair info (pos1, dist_to_pos2, nt1,nt2) into a single int
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
int Bloocoo::decode_index(unsigned int idx, int * pos1, int * dist, int * nt1, int * nt2)
{
    *nt2 = idx & 3;
    idx >>= 2;
    *nt1 = idx & 3;
    idx >>= 2;
    *dist = idx & 15;
    idx >>= 4;
    *pos1 = idx;
    
    return *pos1;
}

