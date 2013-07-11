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

#include <sstream> //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! for test only

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
const char* Bloocoo::STR_NB_ITER_PER_READ          = "-nb-iter";
const char* Bloocoo::STR_NB_VALIDATED_KMERS         = "-nkmer-checked";





/********************************************************************************/
//fonctions for correction
/********************************************************************************/




#define PRINT_DEBUG 0














/********************************************************************************/
//functor for read correction this is the code to correct a single read
/********************************************************************************/
class CorrectReads : public IteratorFunctor
{
public:
    void operator() ( Sequence& s) //no const otherwise error with tKmer.setData
    {
        //if(_bloocoo._seq_num > 116){
        //	return;
        //}
        
        // _bloom   is the bloom filter containing the solid kmers
        // _bloom.contains(kmer)  to ask if it contains a kmer
        
        _bloocoo._corrected_pos.clear();
        
        if(PRINT_DEBUG){	_bloocoo.__badReadStack = "\n\n\n";	}
        
        //char bin2NT[4] = {'A','C','T','G'};
        //char bin2NTrev[4] = {'T','G','A','C'};
        //char binrev[4]    = {2,3,0,1};
        
        //printf("---- new read ----\n");
        
        char * readseq = s.getDataBuffer(); // the nucleotide sequence of the read
        size_t sizeKmer = _bloocoo._kmerSize;
        
        //int nbPasses = _bloocoo._nb_passes_per_read;
        int max_nb_kmers_checked = _bloocoo._nb_kmers_checked;
        
        int nb_checked;
        KmerModel model (sizeKmer);
        KmerModel::Iterator itKmer (model);
        
        
        //int pass;
        //multiple passes per read
        bool continue_correction = true;
        //int nb_errors_cor = 0;
        bool first_gap = true;
        int ii=0;
        int readlen = s.getDataSize();
        
        kmer_type current_kmer;
        kmer_type current_kmer_min;
        //kmer_type previous_kmer = 0;
        kmer_type kmer_begin = 0 ;
        kmer_type kmer_end = 0;
            
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
            										
            
            
            
			if(PRINT_DEBUG){ _bloocoo.print_read_correction_state(&model, s, _bloocoo._seq_num); }
		      
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
                        kmer_end = *kmers[ii-1];
                        
                        if (untrusted_zone_size>1){
                        
                        	
                        	/*
		                    //two-sided conservative correction
		                    // this should be an isolated error, middle of the read
		                    // if error at pos sizeKmer, could work in theory but would need kmer_begin+1
		                    if((untrusted_zone_size == sizeKmer)  &&  (ii != sizeKmer)){
		                        
		                        _bloocoo.__badReadStack += "\t\tTwo sided (hole size k)\n";

		                        nb_errors_cor = twoSidedCorrection(& _bloom , ii-1, readseq,  kmer_begin, kmer_end);
		                        //update_nb_errors_corrected(nb_errors_cor, &_local_nb_errors_corrected, &continue_correction);

		                    }*/
                        
                        	if(nb_errors_cor == 0){
                        	
		                        nb_checked = min(max_nb_kmers_checked, untrusted_zone_size);
		                        nb_checked = max(nb_checked, 1);
		                        
		                        //if first gap, cannot correct from the left side of the gap
		                        if(!first_gap)
		                        {
		                        
		                        	if(PRINT_DEBUG){ _bloocoo.__badReadStack += "\t\tAggressive correction (big hole 1 right)\n"; }
		                            nb_errors_cor = _bloocoo.aggressiveCorrection(ii- untrusted_zone_size -1 , readseq,  kmer_begin, nb_checked , Bloocoo::RIGHT);
		                            
		                        }
		                        
		                        
		                        if(nb_errors_cor == 0){
		                        	if(PRINT_DEBUG){ _bloocoo.__badReadStack += "\t\tAggressive correction (big hole 2 left)\n"; }
		                        	nb_errors_cor = _bloocoo.aggressiveCorrection(ii-1 , readseq,  kmer_end, nb_checked , Bloocoo::LEFT);
		                        }
		                        
				                if(nb_errors_cor == 0){
				                	if(PRINT_DEBUG){ _bloocoo.__badReadStack += "\t\tVote correction (big hole 2)\n"; }
				                	nb_errors_cor = _bloocoo.voteCorrectionInUntrustedZone(ii- untrusted_zone_size, ii-2, readseq, kmers, nb_checked);
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
                    
                    if(trusted_zone_size > 1) // begin of not indexed zone
                    {
                    	kmer_begin = *kmers[ii-1];
                        
                    }
                    untrusted_zone_size ++;
                    trusted_zone_size =0;
                    
                    //end of the read, we should treat this gap here correc snp with trad method
                    if(ii == (readlen-sizeKmer)){
                        nb_checked = min(max_nb_kmers_checked, (int)untrusted_zone_size);
                        if(PRINT_DEBUG){ _bloocoo.__badReadStack += "\t\tAggressive correction (end)\n"; }
                        nb_errors_cor = _bloocoo.aggressiveCorrection(readlen - untrusted_zone_size - sizeKmer , readseq,  kmer_begin, nb_checked , Bloocoo::RIGHT);
                            
                        if(nb_errors_cor == 0){
                        	if(PRINT_DEBUG){ _bloocoo.__badReadStack += "\t\tVote correction (end)\n"; }
							nb_errors_cor = _bloocoo.voteCorrectionInUntrustedZone(readlen - (untrusted_zone_size-1) - sizeKmer , readlen-sizeKmer, readseq, kmers, nb_checked);
                    	}
                        
                    }
                    
                }
                
                _bloocoo.update_nb_errors_corrected(nb_errors_cor, &_local_nb_errors_corrected, &continue_correction);
                
            } // end of kmers iteration over the read
            
            
        }
        
        if(PRINT_DEBUG){ 
        	//printf("%s", _bloocoo.__badReadStack.c_str());
        	_bloocoo.print_read_if_not_fully_corrected(&model, s);
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

    /** We add options specific to this tool. */
    getParser()->add (new OptionOneParam (DSK::STR_KMER_SIZE,               "size of a kmer",   true));
    getParser()->add (new OptionOneParam (DSK::STR_URI_DATABASE,            "database",         true));   // not useful ?
    getParser()->add (new OptionOneParam (DSK::STR_URI_SOLID_KMERS,         "solid kmers file", false));
    getParser()->add (new OptionOneParam (Bloocoo::STR_NB_ITER_PER_READ,    "number of iterations per read", false,"3"));
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
    _nb_passes_per_read = getInput()->getInt (STR_NB_ITER_PER_READ);
    
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
int Bloocoo::twoSidedCorrection(int pos, char *readseq, kmer_type kmer_begin,kmer_type kmer_end)
{
    
    int original_nt = (readseq[pos-1]>>1)&3;

    
    char bin2NT[4] = {'A','C','T','G'};
    char binrev[4] = {2,3,0,1};
    
    KmerModel model (_kmerSize);
    
    //kmer_begin is the last indexed kmer, test its right extension to get the snp
    int good_nt;
    
    //kmer_begin.printASCII(_kmerSize);
    //kmer_end.printASCII(_kmerSize);
    kmer_type temp_kmer, kmer_query;
    
    int nb_alternative = 0;

    for(int nt=0; nt<4; nt++){
		if(nt == original_nt){
			continue;
		}
		
		//Leftmost kmer check
		temp_kmer = codeSeedBin(&model, &kmer_begin, nt, RIGHT);
		kmer_query = min(temp_kmer, revcomp(temp_kmer, _kmerSize));
        //temp_kmer = model.codeSeedRight (kmer_begin, nt, Data::INTEGER);
        //printf("Left most Kmer: "); temp_kmer.printASCII(_kmerSize);
        //temp_kmer.printASCII(_kmerSize);
        
		//kmer is indexed
        if(_bloom->contains(kmer_query)){
            
            //Rightmost kmer check
            //kmer_end_rev =revcomp(kmer_end, _kmerSize);
            //temp_kmer = model.codeSeedRight (kmer_end_rev, binrev[nt], Data::INTEGER);
            temp_kmer = codeSeedBin(&model, &kmer_end, nt, LEFT);
            kmer_query = min(temp_kmer, revcomp(temp_kmer, _kmerSize));
            //printf("Right most Kmer: "); temp_kmer.printASCII(_kmerSize);
            
            if(_bloom->contains(kmer_query)){
                nb_alternative += 1;
                good_nt = nt;
            }
        }
    }
    
    if(nb_alternative == 1)
    {
        return apply_correction(readseq, pos-1, good_nt);
    }
    
    if(PRINT_DEBUG){ __badReadStack += "\t\tfailed\n"; }
    //printf("\t\tfailed\n");
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
int Bloocoo::aggressiveCorrection(int pos, char *readseq, kmer_type kmer_begin,int nb_kmer_check,Direction direction)
{

    //printf("Start agresseive correction\n");
	int votes [4] = {0, 0, 0, 0};

    int good_nt;
	char nt_temp;
    
    char bin2NT[4] = {'A','C','T','G'};
    //char binrev[4]    = {2,3,0,1};
    
    //char revASCII [256];
    //revASCII['A']= 'T';
    //revASCII['T']= 'A';
    //revASCII['C']= 'G';
    //revASCII['G']= 'C';
    
    KmerModel model (_kmerSize);
    
    kmer_type temp_kmer, kmer_query;
    
    
    int original_nt;
    if(direction ==RIGHT)
    {
        original_nt = (readseq[pos+_kmerSize]>>1)&3;
    }
    else
    {
    	original_nt = (readseq[pos-1]>>1)&3;
    }
        
    //printf("original nt: %c", bin2NT[original_nt]);
    for(int nt=0; nt<4; nt++)
    {
		if(nt == original_nt){
			continue;
		}
				
		//printf("\tTest with nt: %c\n", bin2NT[nt]);
		//printf("\t\tkmer begin is: ");
		//kmer_begin.printASCII(_kmerSize);

		temp_kmer = codeSeedBin(&model, &kmer_begin, nt, direction);
        
        //temp_kmer.printASCII(_kmerSize);
        kmer_query = min(temp_kmer, revcomp(temp_kmer, _kmerSize));
        
        if(_bloom->contains(kmer_query)) //kmer is indexed
        {
        	//printf("\t\t\t%c first in bloom  ", bin2NT[nt]);
			votes[nt] += 1;
        }
        else{
        	continue;
        }

				//printf("\t\t");
				//temp_kmer.printASCII(_kmerSize);

		for (int ii=0; ii<nb_kmer_check; ii++)
		{
		    
		    //   printf("checking kmer %i \n",ii);
		    //kmer_begin = temp_kmer;
		    
		    // kmer_begin.printASCII(_kmerSize);
		    
		    if(direction == RIGHT)
		    {
		        nt_temp = readseq[pos+_kmerSize+1+ii];
		    }
		    else
		    {
		        nt_temp = readseq[pos-2-ii];
		    }
		    
		    temp_kmer = codeSeedNT(&model, &temp_kmer, nt_temp, direction);
			
				//printf("\t\t");
				//temp_kmer.printASCII(_kmerSize);
		    
		    kmer_query = min(temp_kmer, revcomp(temp_kmer, _kmerSize));
		    if(_bloom->contains(kmer_query)) //kmer is indexed
		    {
		    	//printf("\t\t\t%c in bloom  ", bin2NT[nt]); temp_kmer.printASCII(_kmerSize);
		        votes[nt] += 1;
				//print_agressive_votes(votes);
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

	int t = (nb_kmer_check*25)/100;
	t = max(1, t);
	if(max_score < t){
		if(PRINT_DEBUG){ __badReadStack += "\t\tfailed\n"; }
		//printf("\t\tfailed\n");
		return 0;
	}
	
	//printf("max: %i\n", max_score);

	//We have the max score in votes, now we have to check if this score
	//is uniq or if more than one nt has the max score
	bool nb_max_score = 0;

	for(int i=0; i<4; i++){
 		if(votes[i] == max_score){
			nb_max_score += 1;
			good_nt = i;
		}
	}

	//printf("good_nt: %i, nb_max_score: %i\n", good_nt, nb_max_score);

	//if max score is uniq we can correct the sequence
	//else if nb_nt_max > 1 then there are more than one candidate for the correction
	//so we cannot determine a good correction, the correction is cancelled.
	if(nb_max_score == 1){

        int pos_corrigee;
        
        if(direction ==RIGHT)
        {
            pos_corrigee= pos+_kmerSize;
        }
        else
        {
            pos_corrigee = pos-1;
            
        }
        
        return apply_correction(readseq, pos_corrigee, good_nt);
	}
    
    if(PRINT_DEBUG){ __badReadStack += "\t\tfailed\n"; }
    
    return 0;
    
}




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
	
	
    //printf("Start voting correction\n");

    int good_nt;
    char bin2NT[4] = {'A','C','T','G'};
    KmerModel model (_kmerSize);
    kmer_type kmer_begin = *kmers[start_pos];
    kmer_type current_kmer, kmer_query, mutated_kmer;

	int nb_column = nb_kmer_check + _kmerSize;
	int votes[nb_column][4];
	for (int i = 0; i < nb_column; i++) {
		for (int nt=0; nt < 4; nt++) {
			votes[i][nt] = 0;
		}
	}
	
	//kmer_begin.printASCII(_kmerSize);
	current_kmer = kmer_begin;

	for(int i=0; i < nb_kmer_check; i++){

		// (pos+1) is the first untrusted base
		int read_pos = (start_pos+1) + i;
		//current_kmer.printASCII(_kmerSize);
		//current_kmer = codeSeedNT(&model, &current_kmer, readseq[read_pos+(_kmerSize-1)], RIGHT);
		current_kmer = model.codeSeedRight(current_kmer, readseq[read_pos+(_kmerSize-1)], Data::ASCII, KMER_DIRECT);
		//current_kmer.printASCII(_kmerSize);
		kmer_query = min(current_kmer, revcomp(current_kmer, _kmerSize));
		if(_bloom->contains(kmer_query)){
		    continue;
		}
		if(is_pos_corrected(read_pos)){
			continue;
		}
				
		for(int kpos=0; kpos < _kmerSize; kpos++){
			int original_nt = (readseq[read_pos+kpos]>>1)&3; //NT2bin ModelAbstract
			//printf("original_nt: %i\n", original_nt);
    		for(int nt=0; nt<4; nt++){
				if(nt == original_nt){
					continue;
				}
				//printf("\t%i\n", nt);
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
				
				kmer_query = min(mutated_kmer, revcomp(mutated_kmer, _kmerSize));
				if(_bloom->contains(kmer_query))
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
		if(PRINT_DEBUG){ __badReadStack += "\t\tfailed\n"; }
		//printf("\t\tfailed\n");
		return 0;
	}
	
	//printf("max vote: %i", maxVote);
	int nb_alternative, pos_corrigee;
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
			pos_corrigee = (start_pos+1) + i;
			/*
			printf("\t\t%i\t%i\t%c\t%c   with   (max_score: %i)  (mini: %i)   (nb_kmer_check: %i)\n",_seq_num,pos_corrigee  ,bin2NT[good_nt], readseq[pos_corrigee], maxVote, t, nb_kmer_check);
			errfile->print("%i\t%i\t%c\t%c\n",_seq_num,pos_corrigee  ,bin2NT[good_nt], readseq[pos_corrigee]);
			readseq[pos_corrigee] = bin2NT[good_nt];*/
			nb_correction += apply_correction(readseq, pos_corrigee, good_nt);
		}
	}

	if(PRINT_DEBUG){
		if(nb_correction == 0){
			__badReadStack += "\t\tfailed because of multiple alternative\n";
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
	if(is_pos_corrected(pos)){
		return 0;
	}
	char bin2NT[4] = {'A','C','T','G'};
	
     if(readseq[pos] =='N') return 0; // dont correct N ?
    
	std::ostringstream oss;
    oss << _seq_num;
    std::string str__seq_num = oss.str();
    std::ostringstream oss2;
    oss2 << pos;
    std::string str_pos = oss2.str();
	if(PRINT_DEBUG){ __badReadStack += "\t\t" + str__seq_num + "\t" + str_pos + "\t" + bin2NT[good_nt] + "\t" + readseq[pos] + "\n"; }
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
	char binrev[4] = {2,3,0,1};
	
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
	char bin2NT[4] = {'A','C','T','G'};
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
	char bin2NT[4] = {'A','C','T','G'};
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
void Bloocoo::print_read_correction_state(KmerModel* model, Sequence& s, int seq_num){
	KmerModel::Iterator itKmer2 (*model);
	kmer_type kmer, kmer_min;
	
	itKmer2.setData (s.getData(),KMER_DIRECT);
	int i3=0;
	//printf("\nread correction state (%i):\n", seq_num);
	//printf("\t%s\n", s.getDataBuffer());
	//printf("\t");
	
	std::ostringstream oss;
    oss << seq_num;


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
bool Bloocoo::is_pos_corrected(int pos){
	return std::find(_corrected_pos.begin(), _corrected_pos.end(), pos) != _corrected_pos.end();
}
