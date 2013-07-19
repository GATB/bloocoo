/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

#include <gatb/system/impl/System.hpp>

#include <gatb/bank/impl/Bank.hpp>
#include <gatb/bank/impl/BankBinary.hpp>
#include <gatb/bank/impl/BankHelpers.hpp>

#include <gatb/kmer/impl/Model.hpp>
#include <gatb/kmer/impl/BloomBuilder.hpp>

#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>
#include <gatb/tools/designpattern/impl/Command.hpp>

#include <gatb/tools/collections/impl/BagFile.hpp>
#include <gatb/tools/collections/impl/BagCache.hpp>
#include <gatb/tools/collections/impl/IteratorFile.hpp>

#include <gatb/tools/misc/impl/Progress.hpp>
#include <gatb/tools/misc/impl/Property.hpp>
#include <gatb/tools/misc/impl/TimeInfo.hpp>

#include <iostream>
#include <map>
#include <math.h>

#include <gatb/tools/math/Integer.hpp>

#include <libgen.h>

#include <Bloocoo.hpp>
#include <DSK.hpp>

// We use the required packages
using namespace std;

using namespace gatb::core::system;
using namespace gatb::core::system::impl;

using namespace gatb::core::bank;
using namespace gatb::core::bank::impl;

using namespace gatb::core::kmer;
using namespace gatb::core::kmer::impl;

using namespace gatb::core::tools::dp;
using namespace gatb::core::tools::dp::impl;

using namespace gatb::core::tools::misc;
using namespace gatb::core::tools::misc::impl;

using namespace gatb::core::tools::collections;
using namespace gatb::core::tools::collections::impl;

using namespace gatb::core::tools::math;

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




#define PRINT_DEBUG 0
#define PRINT_STATS 1













/********************************************************************************/
//functor for read correction this is the code to correct a single read
/********************************************************************************/
class CorrectReads : public IteratorFunctor
{
public:
    void operator() ( Sequence& s) //no const otherwise error with tKmer.setData
    {
        //if(s.getIndex() != 1436){
        //	return;
        //}
        
        //if(s.getIndex() != 4925){
        //	return;
        //}
        // _bloom   is the bloom filter containing the solid kmers
        // _bloom.contains(kmer)  to ask if it contains a kmer
        
        _bloocoo._corrected_pos.clear();
        
        if(PRINT_DEBUG){ _bloocoo.__badReadStack = "\n\n\n"; }
        

        
        //printf("---- new read ----\n");
        
        char * readseq = s.getDataBuffer(); // the nucleotide sequence of the read
        size_t sizeKmer = _bloocoo._kmerSize;
        
        //int nbPasses = _bloocoo._nb_passes_per_read;
        int max_nb_kmers_checked = _bloocoo._nb_kmers_checked;
        
        int nb_checked;
        KmerModel model (sizeKmer);
        KmerModel::Iterator itKmer (model);
        
        
        bool continue_correction = true;
        bool first_gap = true;
        int ii=0;
        int readlen = s.getDataSize();
        
        kmer_type current_kmer;
        kmer_type current_kmer_min;
            
        //for (int pass=0; pass<3; pass++)
        while(continue_correction)
        {
            continue_correction = false;
            first_gap = true;
            ii = 0;
            
            itKmer.setData (s.getData(),KMER_DIRECT);
            
            int untrusted_zone_size = 0;
            int trusted_zone_size = 0;
            
            //Mettre en dehors du while dans une version final (attention dangereux), faire gaffe a ne jamais utiliser un kmer de ce tableau
            //avec un indice > ii
            kmer_type* kmers[readlen-sizeKmer+1]; 
            										
                        
			if(PRINT_DEBUG){ _bloocoo.print_read_correction_state(&model, s); }
		      
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
                
                
                if (!is_last_kmer_indexed_after_hole && _bloom.contains(current_kmer_min)) //kmer is solid
                {
                    
                    trusted_zone_size += 1;
                    
					//beginning of indexed zone
                    if(trusted_zone_size == 2) 
                    {
                       
                        if (untrusted_zone_size>1){
                        
                        	
		                    //two-sided conservative correction
		                    // this should be an isolated error, middle of the read
		                    // if error at pos sizeKmercorrected_pos+1+ii, could work in theory but would need kmer_begin+1
		                    if(((untrusted_zone_size-1) == sizeKmer)  &&  (ii > sizeKmer)){
		                        
                                if(PRINT_DEBUG){ _bloocoo.__badReadStack += "\t\tTwo sided (hole size k)\n";}
		                        nb_errors_cor = _bloocoo.twoSidedCorrection(ii-2, readseq, kmers);

		                    }
                        
                        	if(nb_errors_cor == 0){
                        	
		                        nb_checked = min(max_nb_kmers_checked, untrusted_zone_size-2); // en fait plus utilisé en interne  dapres gaetan ..
		                        nb_checked = max(nb_checked, 0);
		                        
		                        //if first gap, cannot correct from the left side of the gap
		                        //2eme condition: Agressive right peut corriger le debut si le trou est plus long que sizeKmer mais est-ce utile?
		                        if(!first_gap || (first_gap && untrusted_zone_size > sizeKmer))
		                        {
		                        	if(PRINT_DEBUG){ _bloocoo.__badReadStack += "\t\tAggressive correction (big hole 1 right)\n"; }
		                            nb_errors_cor = _bloocoo.aggressiveCorrection(ii-untrusted_zone_size, readseq,  kmers, nb_checked , Bloocoo::RIGHT);
		                        }
		                        
		                        
		                        if(nb_errors_cor == 0){
		                        	if(PRINT_DEBUG){ _bloocoo.__badReadStack += "\t\tAggressive correction (big hole 2 left)\n"; }
		                        	nb_errors_cor = _bloocoo.aggressiveCorrection(ii-2 , readseq,  kmers, nb_checked , Bloocoo::LEFT);
		                        }
		                        
				                if(nb_errors_cor == 0){
				                	if(PRINT_DEBUG){ _bloocoo.__badReadStack += "\t\tVote correction (big hole 2)\n"; }
				                	nb_errors_cor = _bloocoo.voteCorrectionInUntrustedZone(ii- untrusted_zone_size, ii-2, readseq, kmers, nb_checked);
				            	}
                    	
                    			if(nb_errors_cor == 0){
				                	if(PRINT_DEBUG){ _bloocoo.__badReadStack += "\t\tMulti Mutate Vote correction (big hole 2)\n"; }
				                	nb_errors_cor = _bloocoo.multiMutateVoteCorrectionInUntrustedZone(ii- untrusted_zone_size, ii-2, readseq, kmers, nb_checked, 2);
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
                        nb_checked = min(max_nb_kmers_checked, untrusted_zone_size-1);
                        nb_checked = max(nb_checked, 0);
                        
                        if(PRINT_DEBUG){ _bloocoo.__badReadStack += "\t\tAggressive correction (end)\n"; }
                        nb_errors_cor = _bloocoo.aggressiveCorrection(readlen - untrusted_zone_size - sizeKmer + 1, readseq,  kmers, nb_checked , Bloocoo::RIGHT);
                            
                        if(nb_errors_cor == 0){
                        	if(PRINT_DEBUG){ _bloocoo.__badReadStack += "\t\tVote correction (end)\n"; }
							nb_errors_cor = _bloocoo.voteCorrectionInUntrustedZone(readlen - (untrusted_zone_size-1) - sizeKmer , readlen-sizeKmer, readseq, kmers, nb_checked);
                    	}
                        
                        if(nb_errors_cor == 0){
		                	if(PRINT_DEBUG){ _bloocoo.__badReadStack += "\t\tMulti Mutate Vote correction (end)\n"; }
		                	nb_errors_cor = _bloocoo.multiMutateVoteCorrectionInUntrustedZone(readlen - (untrusted_zone_size-1) - sizeKmer, readlen-sizeKmer, readseq, kmers, nb_checked, 2);
		            	}
                    }
                    
                }
                
                _bloocoo.update_nb_errors_corrected(nb_errors_cor, &_local_nb_errors_corrected, &continue_correction);
                
            } // end of kmers iteration over the read
            
            
        }
        
        if(PRINT_DEBUG){ 
        	printf("%s", _bloocoo.__badReadStack.c_str());
        	//_bloocoo.print_read_if_not_fully_corrected(&model, s);
        }
        
        //if (getSynchro()-)
        {
            getSynchro()->lock()  ;
        }
        _outbank.insert (s); //output corrected sequence
        //if (_synchro)  {
            getSynchro()->unlock() ;// }
        
       // _outbank.insert(s); //output corrected sequence
        
        //    printf("%s\n",s.getDataBuffer());
        //    printf("%s\n",readseq);
        
        _bloocoo._seq_num++; //counter 0 based, not thread-safe  //todo : include the sequence number in the sequence type
        
    }

    CorrectReads (Bloom<kmer_type>& bloom, Bank& outbank, Bloocoo & bloocoo, u_int64_t& nb_errors_corrected)
        : _bloom(bloom), _outbank(outbank), _bloocoo(bloocoo),
          _total_nb_errors_corrected (nb_errors_corrected), _local_nb_errors_corrected(0),
          model(_bloocoo._kmerSize), itKmer(model), _synchro(this->newSynchro())
    {}

    ~CorrectReads ()
    {
        /** We increase the global number of corrected errors. */
        __sync_fetch_and_add (&_total_nb_errors_corrected, _local_nb_errors_corrected);
    }

    Bloom<kmer_type>& _bloom; // the bloom containing the solid kmers
    Bank&             _outbank; // the bank cto insert the result : corrected reads
    Bloocoo &         _bloocoo; // the parent bloocoo object

    KmerModel           model;
    KmerModel::Iterator itKmer;

    ISynchronizer* _synchro;
    ISynchronizer* getSynchro ()  { return _synchro; }

    u_int64_t&  _total_nb_errors_corrected;
    u_int64_t   _local_nb_errors_corrected;
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
    _seq_num = 0;

	if(PRINT_STATS){
		__correction_methods_successes[0] = 0; //TwoSided
		__correction_methods_successes[1] = 0; //Agressive
		__correction_methods_successes[2] = 0; //Vote
	}
    	
    /** We add options specific to this tool. */
    getParser()->add (new OptionOneParam (DSK::STR_KMER_SIZE,               "size of a kmer",   true));
    getParser()->add (new OptionOneParam (DSK::STR_URI_DATABASE,            "database",         true));   // not useful ?
    getParser()->add (new OptionOneParam (DSK::STR_URI_SOLID_KMERS,         "solid kmers file", false));
    getParser()->add (new OptionOneParam (Bloocoo::STR_NB_MIN_VALID,    "min number of kmers to valid a correction", false,"2"));
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
    _kmerSize           = getInput()->getInt (DSK::STR_KMER_SIZE);
    _solidFile          = getInput()->getStr (DSK::STR_URI_SOLID_KMERS);
    _nb_kmers_checked   = getInput()->getInt (STR_NB_VALIDATED_KMERS);
    _nb_min_valid = getInput()->getInt (STR_NB_MIN_VALID);
    _max_multimutation_distance = 6;
    
    /*************************************************/
    /** We create a bloom with inserted solid kmers. */
    /*************************************************/
    _bloom = createBloom ();
    LOCAL (_bloom);

    //iterate over initial file
    Bank inbank (getInput()->getStr(DSK::STR_URI_DATABASE));

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
    string prefix = System::file().getBaseName (getInput()->getStr(DSK::STR_URI_DATABASE));

    /** We set the filename as the base name + a specific suffix. */
    string fileName = prefix + string("_corrected.fasta");

    Bank outbank (fileName);

    u_int64_t total_nb_errors_corrected = 0;

    //file with list of errors, for testing purposes
    string ferrfile = prefix + string ("_bloocoo_corr_errs.tab");

    _errfile = System::file().newFile (ferrfile, "wb");

    /*************************************************/
    // We iterate over sequences and correct them
    /*************************************************/
    {
        TIME_INFO (getTimeInfo(), "sequences correction");

     //   CorrectReads fct2(*bloom, outbank, *this,total_nb_errors_corrected) ;
     //   itSeq->iterate (fct2); //
        setDispatcher (new SerialCommandDispatcher());
        
        getDispatcher()->iterate (itSeq,  CorrectReads (*_bloom, outbank, *this, total_nb_errors_corrected)); // not working ?
    }

    /*************************************************/
    // We gather some statistics.
    /*************************************************/
    getInfo()->add (1, "result");
    getInfo()->add (2, "nb errors corrected", "%ld", total_nb_errors_corrected);
    getInfo()->add (2, "corrected file",      fileName);
    if(PRINT_STATS){
    	//getInfo()->add (1, "Correction methods stats:");
    	//getInfo()->add (2, "\tTwoSided:", "%i", __correction_methods_successes[0]);
    	printf("Correction methods stats:\n");
    	printf("\tTwoSided: %i\n", __correction_methods_successes[0]);
    	printf("\tAggressive: %i\n", __correction_methods_successes[1]);
    	printf("\tVote: %i\n", __correction_methods_successes[2]);
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

    u_int64_t solidFileSize = (System::file().getSize(_solidFile) / sizeof (kmer_type));

    u_int64_t estimatedBloomSize = (u_int64_t) ((double)solidFileSize * NBITS_PER_KMER);
    if (estimatedBloomSize ==0 ) { estimatedBloomSize = 1000; }

    /** We create the kmers iterator from the solid file. */
    Iterator<kmer_type>* itKmers = createIterator<kmer_type> (
        new IteratorFile<kmer_type> (_solidFile),
        solidFileSize,
        "fill bloom filter"
    );
    LOCAL (itKmers);

    /** We instantiate the bloom object. */
    BloomBuilder<kmer_type> builder (itKmers, estimatedBloomSize, (int)floorf (0.7*NBITS_PER_KMER), getInput()->getInt(Tool::STR_NB_CORES));
    Bloom<kmer_type>* bloom = builder.build ();

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
int Bloocoo::twoSidedCorrection(int pos, char *readseq, kmer_type* kmers[])
{
	if(!is_pos_correctable(pos, readseq)){
		if(PRINT_DEBUG){ __badReadStack += "\t\t\tfailed (position is already corrected)\n"; }
		return 0;
	}
	
    kmer_type left_most_kmer = *kmers[pos-_kmerSize+1];
    kmer_type right_most_kmer = *kmers[pos];
    
    int original_nt = (readseq[pos]>>1)&3;
    int good_nt;
    kmer_type mutated_kmer, mutated_kmer_min;
    int nb_alternative = 0;

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
    	int nb_correction = apply_correction(readseq, pos, good_nt);
    	if(PRINT_STATS){ __correction_methods_successes[0] += nb_correction; }
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
int Bloocoo::aggressiveCorrection(int pos, char *readseq, kmer_type* kmers[], int nb_kmer_check,Direction direction)
{
	
	kmer_type kmer_begin = *kmers[pos];
    int corrected_pos;
    
    if(direction ==RIGHT){
        corrected_pos = pos+_kmerSize-1;
    }
    else{
        corrected_pos = pos;
    }
        
	if(!is_pos_correctable(corrected_pos, readseq)){
		if(PRINT_DEBUG){ __badReadStack += "\t\t\tfailed (position is already corrected)\n"; }
		return 0;
	}
	
	int votes [4] = {0, 0, 0, 0};

    int good_nt;
	char nt_temp;
    
    
    //char revASCII [256];
    //revASCII['A']= 'T';
    //revASCII['T']= 'A';
    //revASCII['C']= 'G';
    //revASCII['G']= 'C';
    
    KmerModel model (_kmerSize);
    kmer_type current_kmer, current_kmer_min;
    
    int original_nt;
    if(direction ==RIGHT){
        original_nt = (readseq[corrected_pos]>>1)&3;
    }
    else{
    	original_nt = (readseq[corrected_pos]>>1)&3;
    }
        
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
        }
        else{
        	continue;
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
		    
		    current_kmer = codeSeedNT(&model, &current_kmer, nt_temp, direction);
		    
		    current_kmer_min = min(current_kmer, revcomp(current_kmer, _kmerSize));
		    if(_bloom->contains(current_kmer_min)) //kmer is indexed
		    {
		        votes[nt] += 1;
		    }
		    
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
	int t = max(1, nb_kmer_check);
	//int t = 1; //max(1, nb_kmer_check);
	if(max_score < t){
		if(PRINT_DEBUG){ __badReadStack += "\t\t\tfailed (max score %i is too low)\n", max_score; }
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
		int nb_correction = apply_correction(readseq, corrected_pos, good_nt);
    	if(PRINT_STATS){ __correction_methods_successes[1] += nb_correction; }
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
int Bloocoo::voteCorrectionInUntrustedZone(int start_pos, int end_pos, char *readseq, kmer_type* kmers[], int nb_kmer_checked){
	//printf("%i:    %i %i\n", _seq_num, start_pos, end_pos);
	//printf("\n\t");
	//(*kmers[start_pos]).printASCII(_kmerSize);
	//(*kmers[end_pos]).printASCII(_kmerSize);
	
	//start_pos = max(0, start_pos);
	
	int nb_errors_cor = 0;
	
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
		
		nb_errors_cor = voteCorrection(start_pos, readseq, kmers, new_nb_checked);
		start_pos += new_nb_checked;
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
int Bloocoo::voteCorrection(int start_pos, char *readseq, kmer_type* kmers[], int nb_kmer_check)
{
	
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
			if(!is_pos_correctable(read_pos+kpos, readseq)){
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
	int vote, maxVote = 0;
	for (int i = 0; i < nb_column; i++) {
		for (int nt=0; nt < 4; nt++) {
			vote = votes[i][nt];
			if (vote > maxVote) {
				maxVote = vote;
			}
		}
	}

	//int t = (nb_kmer_check*25)/100;
	//t = max(2, t);
	int t = max(1, nb_kmer_check);
	if(maxVote < t){
		if(PRINT_DEBUG){ __badReadStack += "\t\t\tfailed\n"; }
		return 0;
	}
	
	//printf("max vote: %i", maxVote);
	int nb_alternative, corrected_pos;
	int nb_correction = 0;

	
	
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
			int nb_cor = apply_correction(readseq, corrected_pos, good_nt);
			nb_correction += nb_cor;
	    	if(PRINT_STATS){ __correction_methods_successes[2] += nb_cor; }
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
int Bloocoo::apply_correction(char *readseq, int pos, int good_nt){
	if(!is_pos_correctable(pos, readseq)){
		if(PRINT_DEBUG){ __badReadStack += "\t\t\tFailed in apply_correction (non correctable pos)\n"; }
		return 0;
	}
	
    // char bin2NT[4] = {'A','C','T','G'};


	if(PRINT_DEBUG){
		std::ostringstream oss;
		oss << _seq_num;
		std::string str__seq_num = oss.str();
		std::ostringstream oss2;
		oss2 << pos;
		std::string str_pos = oss2.str();
		__badReadStack += "\t\t\t" + str__seq_num + "\t" + str_pos + "\t" + bin2NT[good_nt] + "\t" + readseq[pos] + "\n";
	}
	//printf("\t\t%i\t%i\t%c\t%c\n",_seq_num,pos,bin2NT[good_nt],readseq[pos]);
	
    _errfile->print("%i\t%i\t%c\t%c\n",_seq_num,pos,bin2NT[good_nt],readseq[pos]);
    readseq[pos] = bin2NT[good_nt]; //correc sequence in ram
    // printf("error found pos %i  read %i\n",ii, _bloocoo._seq_num);
    _corrected_pos.push_back(pos);
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
kmer_type Bloocoo::codeSeedBin(KmerModel* model, kmer_type* kmer, int nt, Direction direction){
	kmer_type temp_kmer;
	
	if(direction == RIGHT){
		temp_kmer = model->codeSeedRight(*kmer, nt, Data::INTEGER, KMER_DIRECT);
	}
	else{
		kmer_type kmer_end_rev = revcomp(*kmer, _kmerSize);
		temp_kmer = model->codeSeedRight (kmer_end_rev, binrev[nt], Data::INTEGER, KMER_DIRECT);
		temp_kmer = revcomp(temp_kmer, _kmerSize);
	}
	return temp_kmer; //min(temp_kmer, revcomp(temp_kmer, _kmerSize));
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
kmer_type Bloocoo::codeSeedNT(KmerModel* model, kmer_type* kmer, char nt, Direction direction){
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
bool Bloocoo::is_pos_correctable(int pos, char* readseq){
	bool N_at_pos = (readseq[pos] == 'N');
	bool is_pos_corrected = (std::find(_corrected_pos.begin(), _corrected_pos.end(), pos) != _corrected_pos.end());
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
int Bloocoo::multiMutateVoteCorrectionInUntrustedZone(int start_pos, int end_pos, char *readseq, kmer_type* kmers[], int nb_kmer_checked, int max_nb_mutation){
	//printf("%i:    %i %i\n", _seq_num, start_pos, end_pos);
	//printf("\n\t");
	//(*kmers[start_pos]).printASCII(_kmerSize);
	//(*kmers[end_pos]).printASCII(_kmerSize);
	
	//start_pos = max(0, start_pos);
	//printf("\n\tMulti mutate start !!!!!!!!!!\n");
	
	int nb_errors_cor = 0;
	
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
		
		nb_errors_cor = multiMutateVoteCorrection(start_pos, readseq, kmers, new_nb_checked, max_nb_mutation);
		start_pos += new_nb_checked;
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
int Bloocoo::multiMutateVoteCorrection(int start_pos, char *readseq, kmer_type* kmers[], int nb_kmer_check, int max_nb_mutation)
{

	int nb_column = nb_kmer_check + _kmerSize - 1;
	std::map<std::string, int> votes;
	
	for(int i=0; i < nb_kmer_check; i++){
		kmer_type current_kmer = *kmers[start_pos+i];
		multiMutateVoteCorrectionRec(start_pos, i, current_kmer, readseq, kmers, nb_kmer_check, 0, max_nb_mutation, 0, "", votes);
	}
	
	//Search max vote in the hash
	int nbMax = 0;
	int vote, maxVote = 0;
	std::string good_mutations;
	
	std::map<std::string, int>::iterator iter;
	
	for (iter = votes.begin(); iter != votes.end(); ++iter) {
		//printf("%s\n", iter->first.c_str());
		int vote = iter->second;
		if (vote > maxVote) {
			maxVote = vote;
		}
	}
	
	int t = max(3, nb_kmer_check-1);
	if(maxVote < t){
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
			nb_correction += nb_cor;
			i = 0;
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
int Bloocoo::multiMutateVoteCorrectionRec(int start_pos, int kmer_offset, kmer_type current_kmer, char *readseq, kmer_type* kmers[], int nb_kmer_check, int kmer_index, int max_nb_mutation, int current_nb_mutation, std::string mutations, std::map<std::string, int>& votes){
	
	int read_pos = start_pos + kmer_offset;
	
	//kmer_type current_kmer_min = min(current_kmer, revcomp(current_kmer, _kmerSize));
	//if(_bloom->contains(current_kmer_min)){
	//    return 0;
	//}
	//if(!is_pos_correctable(read_pos, readseq)){
	//	continue;
	//}
	
	
	//Determine the distance between 2 mutations
	int end_kpos = _kmerSize;
	if(current_nb_mutation == max_nb_mutation-1){
		end_kpos = std::min(kmer_index+_max_multimutation_distance, (int)_kmerSize);
		//printf("%i %i\n", kmer_index, end_kpos);
	}
	
	//Iterating on kmer positions from kmer_index to the max mutate distance
	//If the max number of mutation is not reached, we iterate over all the kmer positions
	for(int kpos=kmer_index; kpos < end_kpos; kpos++){
		
		//Don't vote if the current position is not correctable
		if(!is_pos_correctable(read_pos+kpos, readseq)){
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
			
			//Add this new mutation ("pos nt") to the string that stores the mutations chain. (Ex: "60 A 12 B...")
			//_oss.clear();
			_oss.str("");
			_oss << kmer_offset+kpos;
			std::string mutations_copy = mutations + " " + _oss.str();
			//_oss.clear();
			_oss.str("");
			_oss << nt;
			mutations_copy += " " + _oss.str();
				
			//If the max number of mutations is reached, we can start voting for the mutated kmers
			if(current_nb_mutation == max_nb_mutation-1){
				/*
				printf("kmer index: %i\n", kpos);
				printf("\tcurent kmer:  "); current_kmer.printASCII(_kmerSize);
				printf("\tmutated kmer:  "); mutated_kmer.printASCII(_kmerSize);
				printf("\n");*/
				
				kmer_type mutated_kmer_min = min(mutated_kmer, revcomp(mutated_kmer, _kmerSize));
				if(_bloom->contains(mutated_kmer_min)){
					votes[mutations_copy] += 1;
				}
			}
			//If the max number of mutations is not reached, we have to recurcivelly continue to applied mutations before voting
			//The kmer_index and the current number of mutation is incremented 
			else{
				multiMutateVoteCorrectionRec(start_pos, kmer_offset, mutated_kmer, readseq, kmers, nb_kmer_check, kpos+kmer_index+1, max_nb_mutation, current_nb_mutation+1, mutations_copy, votes);
			}
		}
	}
}







