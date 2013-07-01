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
const char* Bloocoo::STR_NB_ITER_PER_READ          = "-nb-iter";
const char* Bloocoo::STR_NB_VALIDATED_KMERS         = "-nkmer-checked";





/********************************************************************************/
//fonctions for correction
/********************************************************************************/
int voteCorrection(Bloom<kmer_type> * _bloom , int pos, char *readseq, IFile* errfile, kmer_type kmer_begin, size_t sizeKmer, int numseq, int nb_kmer_check);

typedef enum direction
{
    LEFT,
    RIGHT
} direction_t;


#define NDEB -1

//// print_agressive_votes
//Debug function for agressive correction, print the votes tab
void print_agressive_votes(int votes[4]){
	char bin2NT[4] = {'A','C','T','G'};
	printf("#############################################\n");
	for(int i=0; i<4; i++){
		printf("%c: %i\n", bin2NT[i], votes[i]);
	}
	printf("#############################################\n");
}

//print_votes
//Debug function for vote correction, print the votes matrix
void print_votes(int votes[][4], int nb_column){
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

//// mutate_kmer
// fonction to mutate kmer : takes kmer, pos (de la fin ,0-based), and nt (nt = 0,1,2ou 3) 
// par exemple  kmer ,2 , C    avec kmer =  AAAAAAAAAA
//
// return : 
// AAAAAAACAA
void mutate_kmer(kmer_type * kmer, int pos, char nt)
{
    kmer_type reset_mask =  ~((kmer_type)3 << (pos*2));
    kmer_type set_mask =  ((kmer_type)nt) << (pos*2);
    *kmer = (*kmer & reset_mask )  | set_mask;
}

//// twoSidedCorrection
// Method 1 for read correction
// Correct an error in an untrusted zone of size=kmerSize
int twoSidedCorrection(Bloom<kmer_type> * _bloom , int pos, char *readseq, IFile* errfile, kmer_type kmer_begin,kmer_type kmer_end, size_t  sizeKmer, int numseq)
{
    
    int original_nt = (readseq[pos-1]>>1)&3;

    
    char bin2NT[4] = {'A','C','T','G'};
    char binrev[4]    = {2,3,0,1};
    
    KmerModel model (sizeKmer);
    
    //kmer_begin is the last indexed kmer, test its right extension to get the snp
    int good_nt;
    
    //kmer_begin.printASCII(sizeKmer);
    //kmer_end.printASCII(sizeKmer);
    kmer_type temp_kmer, kmer_end_rev;
    
    int nb_alternative = 0;

    for(int nt=0; nt<4; nt++)
    {
		if(nt == original_nt){
			continue;
		}
		
        temp_kmer = model.codeSeedRight (kmer_begin, nt, Data::INTEGER);
        
        //  temp_kmer.printASCII(sizeKmer);
        
        if(_bloom->contains(temp_kmer)) //kmer is indexed
        {
            //ref is the last nt of the first non indexed kmer on the reference genome
            //contig, pos , ref , snp
            //  fprintf(snp_file,"%i  %i  %c %c  \n",j,i-1,bin2NT[first_unindexed & 3],bin2NT[nt]);
            
            //rightmost kmer check
            kmer_end_rev =revcomp(kmer_end, sizeKmer);
            temp_kmer = model.codeSeedRight (kmer_end_rev, binrev[nt], Data::INTEGER);
            
            if( _bloom->contains(temp_kmer))
            {
                nb_alternative += 1;
                good_nt = nt;
            }
        }
    }
    
    if(nb_alternative == 1)
    {
    	printf("\t\t%i\t%i\t%c\t%c\n",numseq,pos-1,bin2NT[good_nt],readseq[pos-1]);
    	
        errfile->print("%i\t%i\t%c\t%c\n",numseq,pos-1,bin2NT[good_nt],readseq[pos-1]);
        readseq[pos-1] = bin2NT[good_nt]; //correc sequence in ram
        // printf("error found pos %i  read %i\n",ii, _bloocoo._seq_num);
        
        return 1;
    }
    
    printf("\t\tfailed\n");
    return 0;
    
}



//// aggressiveCorrection
// Method 2 for read correction
// This method can correct an error in an untrusted zone of size > kmerSize (close errors)
//
// nb_kmer_check the number of additional kmers checked
// pos = position of the solid kmer next to the untrusted zone
// kmer_begin  = solid kmer at pos
int aggressiveCorrection(Bloom<kmer_type> * _bloom , int pos, char *readseq, IFile* errfile, kmer_type kmer_begin, size_t  sizeKmer, int numseq,int nb_kmer_check,direction_t orientation )
{

    //printf("Start agresseive correction\n");
	int votes [4] = {0, 0, 0, 0};

    int good_nt;
	char nt_temp;
    
    char bin2NT[4] = {'A','C','T','G'};
    char binrev[4]    = {2,3,0,1};
    
    char revASCII [256];
    revASCII['A']= 'T';
    revASCII['T']= 'A';
    revASCII['C']= 'G';
    revASCII['G']= 'C';
    
    KmerModel model (sizeKmer);
    
    
    kmer_type temp_kmer,kmer_end_rev, kmer_query;
    
    
    int original_nt;
    if(orientation ==RIGHT)
    {
        original_nt = (readseq[pos+sizeKmer]>>1)&3;
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
		//kmer_begin.printASCII(sizeKmer);

        if(orientation == RIGHT)
        {
            temp_kmer = model.codeSeedRight (kmer_begin, nt, Data::INTEGER,KMER_DIRECT);
        }
        else
        {
            kmer_end_rev =revcomp(kmer_begin, sizeKmer);
            temp_kmer = model.codeSeedRight (kmer_end_rev, binrev[nt], Data::INTEGER,KMER_DIRECT);
            temp_kmer =revcomp(temp_kmer, sizeKmer);
        }

        
        kmer_query = min(temp_kmer, revcomp(temp_kmer, sizeKmer));
        
        if(_bloom->contains(kmer_query)) //kmer is indexed
        {
			votes[nt] += 1;
        }

				//printf("\t\t");
				//temp_kmer.printASCII(sizeKmer);

		for (int ii=0; ii<nb_kmer_check; ii++)
		{
		    
		    //   printf("checking kmer %i \n",ii);
		    //kmer_begin = temp_kmer;
		    
		    // kmer_begin.printASCII(sizeKmer);
		    
		    if(orientation == RIGHT)
		    {
		        nt_temp = readseq[pos+sizeKmer+1+ii];
		    }
		    else
		    {
		        nt_temp = readseq[pos-2-ii];
		    }
		    
		    if(orientation == RIGHT)
		    {
		        temp_kmer = model.codeSeedRight (temp_kmer, nt_temp, Data::ASCII,KMER_DIRECT);
		    }
		    else
		    {
		        kmer_end_rev =revcomp(temp_kmer, sizeKmer);
		        temp_kmer = model.codeSeedRight (kmer_end_rev, revASCII[nt_temp], Data::ASCII,KMER_DIRECT);
		        
		        temp_kmer =revcomp(temp_kmer, sizeKmer);
		        
		    }
			
				//printf("\t\t");
				//temp_kmer.printASCII(sizeKmer);
		    
		    kmer_query = min(temp_kmer, revcomp(temp_kmer, sizeKmer));
		    if(_bloom->contains(kmer_query)) //kmer is indexed
		    {
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

	int t = (nb_kmer_check*75)/100;
	t = max(1, t);
	if(max_score < t){
		printf("\t\tfailed\n");
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
        
        if(orientation ==RIGHT)
        {
            pos_corrigee= pos+sizeKmer;
        }
        else
        {
            pos_corrigee = pos-1;
            
        }
        
		//if(bin2NT[good_nt] == readseq[pos_corrigee]){ // a supprimer !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		//	return 0;
		//}
		
        if(numseq ==NDEB) printf("Corr validee pos %i good nt %c \n",pos_corrigee,bin2NT[good_nt]);
        //printf("Corr validee pos %i good nt %c \n",pos_corrigee,bin2NT[good_nt]);
        errfile->print("%i\t%i\t%c\t%c\n",numseq,pos_corrigee  ,bin2NT[good_nt], readseq[pos_corrigee]);
        
        printf("\t\t%i\t%i\t%c\t%c   with   (max_score: %i)  (mini: %i)   (nb_kmer_check: %i)\n",numseq,pos_corrigee  ,bin2NT[good_nt], readseq[pos_corrigee], max_score, t, nb_kmer_check);
        
        readseq[pos_corrigee ]=bin2NT[good_nt]; //correc sequence in ram
        
        return 1;
	}
    
    printf("\t\tfailed\n");
    return 0;
    
}


//// voteCorrection
// Method 3 for read correction
// This method is used if the method 2 fails to correct an error
//
// nb_kmer_check the number of additional kmers checked
// pos: the last trusted position in readseq before the untrusted zone
// kmer_begin: solid kmer at pos, it's the last trusted kmer before the untrusted zone
int voteCorrection(Bloom<kmer_type> * _bloom , int pos, char *readseq, IFile* errfile, kmer_type kmer_begin, size_t sizeKmer, int numseq, int nb_kmer_check)
{
	
    //printf("Start voting correction\n");

    int good_nt;
    char bin2NT[4] = {'A','C','T','G'};
    KmerModel model (sizeKmer);
    kmer_type current_kmer, kmer_query, muted_kmer;

	int nb_column = nb_kmer_check + (sizeKmer - 1);
	int votes[nb_column][4];
	for (int i = 0; i < nb_column; i++) {
		for (int nt=0; nt < 4; nt++) {
			votes[i][nt] = 0;
		}
	}
	
	//kmer_begin.printASCII(sizeKmer);
	current_kmer = kmer_begin;

	for(int i=0; i < nb_kmer_check; i++){

		// (pos+1) is the first untrusted base
		int read_pos = (pos+1) + i;
		current_kmer = model.codeSeedRight(current_kmer, readseq[read_pos+(sizeKmer-1)], Data::ASCII, KMER_DIRECT);
		
		//kmer_query = min(current_kmer, revcomp(current_kmer, sizeKmer));
		if(_bloom->contains(current_kmer))
		{
		    continue;
		}
				
		//current_kmer.printASCII(sizeKmer);

		for(int kpos=0; kpos < sizeKmer; kpos++){
			int original_nt = (readseq[read_pos+kpos]>>1)&3; //NT2int ModelAbstract
			//printf("original_nt: %i\n", original_nt);
    		for(int nt=0; nt<4; nt++){
				if(nt == original_nt){
					continue;
				}
				//printf("\t%i\n", nt);
				muted_kmer = current_kmer;
				mutate_kmer(&muted_kmer, sizeKmer-1-kpos, nt);
				//muted_kmer.printASCII(sizeKmer);
				
				/*
				printf("mutate at %i with %c\n", sizeKmer-1-kpos, bin2NT[nt]);
				printf("curre_kmer: ");
				current_kmer.printASCII(sizeKmer);
				printf("muted_kmer: ");
				muted_kmer.printASCII(sizeKmer);
				printf("\n");
				*/
				
				kmer_query = min(muted_kmer, revcomp(muted_kmer, sizeKmer));
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

	int t = (nb_kmer_check*75)/100;
	t = max(1, t);
	if(maxVote < t){
		printf("\t\tfailed\n");
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
			pos_corrigee = (pos+1) + i;
			printf("\t\t%i\t%i\t%c\t%c   with   (max_score: %i)  (mini: %i)   (nb_kmer_check: %i)\n",numseq,pos_corrigee  ,bin2NT[good_nt], readseq[pos_corrigee], maxVote, t, nb_kmer_check);
			errfile->print("%i\t%i\t%c\t%c\n",numseq,pos_corrigee  ,bin2NT[good_nt], readseq[pos_corrigee]);
			readseq[pos_corrigee] = bin2NT[good_nt];
			nb_correction += 1;
		}
	}

    return nb_correction;
    
}

void update_nb_errors_corrected(int nb_errors_corrected, u_int64_t* _local_nb_errors_corrected, bool* continue_correction){
    *_local_nb_errors_corrected += nb_errors_corrected;
    if(nb_errors_corrected > 0){
    	*continue_correction = true;
    }
}




/********************************************************************************/
//functor for read correction this is the code to correct a single read
/********************************************************************************/
class CorrectReads : public IteratorFunctor
{
public:
    void operator() ( Sequence& s) //no const otherwise error with tKmer.setData
    {
        
        // _bloom   is the bloom filter containing the solid kmers
        // _bloom.contains(kmer)  to ask if it contains a kmer
        
        char bin2NT[4] = {'A','C','T','G'};
        char bin2NTrev[4] = {'T','G','A','C'};
        char binrev[4]    = {2,3,0,1};
        
        //printf("---- new read ----\n");
        
        char * readseq = s.getDataBuffer (); // the nucleotide sequence of the read
        size_t   sizeKmer = _bloocoo._kmerSize;
        IFile*  errfile = _bloocoo._errfile;
        
        int nbPasses = _bloocoo._nb_passes_per_read;
        int max_nb_kmers_checked = _bloocoo._nb_kmers_checked;
        
        int nb_checked;
        KmerModel model (sizeKmer);
        KmerModel::Iterator itKmer (model);
        
        //int pass;
        //multiple passes per read
        bool continue_correction = true;
        
        //for (int pass=0; pass<3; pass++)
        while(continue_correction)
        {
            continue_correction = false;
            //if(_bloocoo._seq_num ==NDEB) printf("---pass  %i---\n",pass);
            
            kmer_type kmer_end_rev,temp_kmer;
            
            kmer_type graine, graine_revcomp;
            kmer_type current_kmer ;
            kmer_type previous_kmer = 0;
            kmer_type kmer_begin = 0 ;
            kmer_type first_unindexed = 0;
            kmer_type kmer_end = 0;
            bool first_gap = true;
            int nb_errors_cor = 0;
            //int nb_alternatives_found = 0;
            kmer_type kmer_test, kmer_test_min;
            //for each kmer in this Sequence
            
            // sets itKmer to iterate over all the kmers of the read s
            itKmer.setData (s.getData(),KMER_DIRECT);
            uint64_t tai_not_indexed =0;
            uint64_t tai_previous_break =0;
            uint64_t tai_indexed = 0;
            int readlen = s.getDataSize();
            //bool check = false;
            
            int ii=0;
            
            /////////////////////////////////////////////// print etat read
        	KmerModel::Iterator itKmer2 (model);
        	itKmer2.setData (s.getData(),KMER_DIRECT);
        	int i3=0;
	        printf("\nread state (%i):\n", _bloocoo._seq_num);
	        //printf("\t");
	        //for(int i=0; i<readlen; i++){
	        //	printf("%i", i%10);
	        //}
	        printf("\t%s\n", readseq);
	        printf("\t");
	        for (itKmer2.first(); !itKmer2.isDone(); itKmer2.next(),i3++){
	        	kmer_test = *itKmer2;
            	kmer_test_min = std::min(  revcomp (kmer_test, sizeKmer),kmer_test );
            	
	        	if (_bloom.contains(kmer_test_min)){
	        		printf("1");
	        	}
	        	else{
	        		printf("0");
	        	}
	        }
	        printf("\n");
			///////////////////////////////////////////////
		      
            // We iterate the kmers of this sequence
            for (itKmer.first(); !itKmer.isDone(); itKmer.next(),ii++)
            {
            
                //graine est kmer qui commenece à ii,  dans le sens du read
                //current_kmer est le min de graine avec son revcomp
                
                graine = *itKmer;
                current_kmer = std::min(  revcomp (graine, sizeKmer),graine );
                
                
                if (_bloom.contains(current_kmer)) //kmer is solid
                {
                    tai_indexed++;
                    
                    
                    if(tai_indexed==1) //beginning of indexed zone
                    {
                        kmer_end = graine; // kmer_end should be first kmer indexed after a hole
                        
                        //printf("taille not indexed: %i\n", tai_not_indexed);
                        
                        //two-sided conservative correction
                        // this should be an isolated error, middle of the read
                        // if error at pos sizeKmer, could work in theory but would need kmer_begin+1
                        if((tai_not_indexed == sizeKmer)  &&  (ii != sizeKmer)){
                            if(_bloocoo._seq_num ==NDEB){printf(".. two sided  pos %i \n",ii);kmer_begin.printASCII(sizeKmer);}
                            printf("\t\ttwo sided %i\n", ii);
                            nb_errors_cor = twoSidedCorrection(& _bloom , ii, readseq,errfile,  kmer_begin, kmer_end,   sizeKmer,  _bloocoo._seq_num);
                            update_nb_errors_corrected(nb_errors_cor, &_local_nb_errors_corrected, &continue_correction);
                            if(nb_errors_cor > 0){
                            	break;
                            }
                            
                        }
                        // an error at the beginning of the read : tai_not_indexed is smaller even for a single snp
                        //and correct it from the right
                        else if ((tai_not_indexed < sizeKmer) && first_gap && tai_not_indexed>0){
                            
                            printf("\t\tagressive begin %i\n", ii);
                            nb_checked = min(max_nb_kmers_checked, (int)tai_not_indexed-1);
                            nb_errors_cor = aggressiveCorrection(& _bloom , ii , readseq, errfile,  kmer_end, sizeKmer, _bloocoo._seq_num, nb_checked ,LEFT);
                            update_nb_errors_corrected(nb_errors_cor, &_local_nb_errors_corrected, &continue_correction);
                            
                            if(nb_errors_cor == 0){
                            	//printf("\t\tvote correction begin %i\n", ii);
                            	//nb_errors_cor = voteCorrection(& _bloom , ii, readseq, errfile,  kmer_end, sizeKmer, _bloocoo._seq_num, nb_checked);
                            	//update_nb_errors_corrected(nb_errors_cor, &_local_nb_errors_corrected, &continue_correction);
                            }
                            
                            if(nb_errors_cor > 0){
                            	break;
                            }
                            
                            if(_bloocoo._seq_num ==NDEB){printf(".. aggressive correc begin  pos %i \n",ii);kmer_end.printASCII(sizeKmer);}
                            
                        }
                        // use aggressive correct for close errors or error at beginning
                        else if (tai_not_indexed > sizeKmer){ //|| (sizeKmer==ii && (tai_not_indexed == sizeKmer))){
                            nb_checked = min(max_nb_kmers_checked, (int) (tai_not_indexed-sizeKmer -1));
                            
                            if(_bloocoo._seq_num ==NDEB)
                            {printf(".. aggressive correc pos %i left right  %lli %li \n",ii,tai_not_indexed,sizeKmer);kmer_end.printASCII(sizeKmer);}
                            
                            nb_errors_cor = 0;
                            
                            //if first gap, cannot correct from the left side of the gap
                            if(!first_gap)
                            {
                            	printf("\t\tagressive big hole 1 %i\n", ii);
                                nb_errors_cor = aggressiveCorrection(& _bloom , ii- tai_not_indexed -1 , readseq, errfile,  kmer_begin, sizeKmer, _bloocoo._seq_num, nb_checked ,RIGHT);
                            	update_nb_errors_corrected(nb_errors_cor, &_local_nb_errors_corrected, &continue_correction);
                            
                            	if(nb_errors_cor == 0){
                            		//printf("\t\tvote correction big hole 1 %i\n", ii);
                            		//nb_errors_cor = voteCorrection(& _bloom , ii- tai_not_indexed -1, readseq, errfile,  kmer_begin, sizeKmer, _bloocoo._seq_num, nb_checked);
                            		//update_nb_errors_corrected(nb_errors_cor, &_local_nb_errors_corrected, &continue_correction);
                            	}
                            	
                            }
                            
                            if(nb_errors_cor == 0){
                            	printf("\t\tagressive big hole 2 %i\n", ii);
                            	nb_errors_cor = aggressiveCorrection(& _bloom , ii , readseq, errfile,  kmer_end, sizeKmer, _bloocoo._seq_num, nb_checked ,LEFT);
                            	update_nb_errors_corrected(nb_errors_cor, &_local_nb_errors_corrected, &continue_correction);
                            }
                            
                        	if(nb_errors_cor == 0){
                        		//printf("\t\tvote correction big hole 2 %i\n", ii);
                        		//nb_errors_cor = voteCorrection(& _bloom , ii, readseq, errfile,  kmer_end, sizeKmer, _bloocoo._seq_num, nb_checked);
                        		//update_nb_errors_corrected(nb_errors_cor, &_local_nb_errors_corrected, &continue_correction);
                        	}
	                        if(nb_errors_cor > 0){
	                        	break;
	                        }
                        }
                        
                    }
                    
                    
                    first_gap = false;
                    
                    if(_bloocoo._seq_num ==NDEB) {printf("1 pos %i  : ",ii);graine.printASCII(sizeKmer); }
                    
                    ////////////traitement d'un faux positif du bloom isolé
                    //   if(tai_indexed==1) tai_previous_break = tai_not_indexed; // begin of indexed zone    ( if tai_not_indexed ?)
                    if (tai_indexed > 1) // do not reset  tai_not_indexed if a single positive (might be a FP)
                    {
                        tai_not_indexed = 0; //reset previous hole size
                    }
                    if (tai_indexed==1)// begin of indexed zone
                    {
                        if(_bloocoo._seq_num ==NDEB) printf("not indexed %lli pos %i\n",tai_not_indexed,ii);
                        tai_not_indexed++; //  treat a single solidkmer  as an erroneous kmer
                    }
                    //////////
                    
                    //version sans traitement du FP du bloom
					//if(tai_indexed==1){
					//	tai_not_indexed = 0;
					//}
                    
                }
                else //kmer contains an error
                {
                    if(_bloocoo._seq_num ==NDEB) {printf("0 pos %i  : ",ii);graine.printASCII(sizeKmer); }
                    
                    if(tai_indexed > 1) // begin of not indexed zone
                    {
                        kmer_begin = previous_kmer ; //kmer_begin is the last kmer indexed before the hole
                        first_unindexed = graine;
                        if(_bloocoo._seq_num ==NDEB) printf("indexed %lli pos %i\n",tai_indexed,ii);
                        
                    }
                    tai_not_indexed ++;
                    tai_indexed =0;
                    
                    //end of the read, we should treat this gap here correc snp with trad method
                    if(ii == (readlen-sizeKmer)){
                        nb_checked = min(max_nb_kmers_checked, (int)tai_not_indexed-1);
                        
                        if(_bloocoo._seq_num ==NDEB) {printf("Correc right extrem %i ",ii);kmer_begin.printASCII(sizeKmer); }
                        printf("\t\tagressive end %i\n", ii);
                        //printf("taille not indexed: %i\n", tai_not_indexed);
                        nb_errors_cor = aggressiveCorrection(& _bloom , readlen - tai_not_indexed  - sizeKmer , readseq, errfile,  kmer_begin, sizeKmer, _bloocoo._seq_num, nb_checked ,RIGHT);
                        update_nb_errors_corrected(nb_errors_cor, &_local_nb_errors_corrected, &continue_correction);
                            
                        if(nb_errors_cor == 0){
                        	//printf("\t\tvote correction %i\n", ii);
                        	//nb_errors_cor = voteCorrection(& _bloom , readlen - tai_not_indexed  - sizeKmer, readseq, errfile,  kmer_begin, sizeKmer, _bloocoo._seq_num, nb_checked);
                    		//update_nb_errors_corrected(nb_errors_cor, &_local_nb_errors_corrected, &continue_correction);
                    	}
                        if(nb_errors_cor > 0){
                        	break;
                        }
                        
                    }
                    
                }
                
                previous_kmer = graine; // should be the kmer in its original sense
            } // end of kmers iteration over the read
            
            if(_bloocoo._seq_num ==NDEB)
            {
                if(tai_indexed)
                    printf("indexed %lli pos %i\n",tai_indexed,ii);
                else
                    printf("not indexed %lli pos %i\n",tai_not_indexed,ii);
            }
            
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
    Bloom<kmer_type>* bloom = createBloom ();
    LOCAL (bloom);

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
        
        getDispatcher()->iterate (itSeq,  CorrectReads (*bloom, outbank, *this, total_nb_errors_corrected)); // not working ?
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
    BloomBuilder builder (itKmers, estimatedBloomSize, (int)floorf (0.7*NBITS_PER_KMER), getInput()->getInt(Tool::STR_NB_CORES));
    Bloom<kmer_type>* bloom = builder.build ();

    /** We return the created bloom filter. */
    return bloom;
}

