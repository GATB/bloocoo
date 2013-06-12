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

#include <omptl/omptl_numeric>
#include <omptl/omptl_algorithm>

#include <omp.h>

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
//functor for read correction this is the code to correct a single read
/********************************************************************************/
struct CorrectReads
{
    void operator() ( Sequence& s) //no const otherwise error with tKmer.setData
    {
#if 1
        // _bloom   is the bloom filter containing the solid kmers
        // _bloom.contains(kmer)  to ask if it contains a kmer

        char bin2NT[4] = {'A','C','T','G'};
        char bin2NTrev[4] = {'T','G','A','C'};
        char binrev[4]    = {2,3,0,1};

        //_bloocoo._seq_num++; //counter 1 based, not thread-safe
        //printf("---- new read ----\n");

        char * readseq = s.getDataBuffer (); // the nucleotide sequence of the read
        size_t   sizeKmer = _bloocoo._kmerSize;

        int pass;
        //multiple passes per read
        for (pass=0; pass<2; pass++)
        {

            kmer_type graine, graine_revcomp;
            kmer_type current_kmer ;
            kmer_type previous_kmer = 0;
            kmer_type kmer_begin = 0 ;
            kmer_type first_unindexed = 0;
            kmer_type kmer_end = 0;
            bool first_gap = true;

            //for each kmer in this Sequence

            // sets itKmer to iterate over all the kmers of the read s
            itKmer.setData (s.getData(),KMER_DIRECT);

            uint64_t tai_not_indexed =0;
            uint64_t tai_previous_break =0;
            uint64_t tai_indexed = 0;
            int readlen = s.getDataSize();
            bool check = false;

            // We iterate the kmers of this sequence
            int ii=0;
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

                        if(tai_not_indexed == (sizeKmer))  // this should be an isolated error, middle of the read
                        {

                            //kmer_begin is the last indexed kmer, test its right extension to get the snp
                            char nt;
                            check = false;

                            kmer_type temp_kmer;
                            // kmer_begin.printASCII(sizeKmer);

                            for(nt=0; nt<4; nt++)
                            {
                                temp_kmer = model.codeSeedRight (kmer_begin, nt, Data::INTEGER);

                                //  temp_kmer.printASCII(sizeKmer);

                                if(_bloom.contains(temp_kmer)) //kmer is indexed
                                {
                                    //ref is the last nt of the first non indexed kmer on the reference genome
                                    //contig, pos , ref , snp
                                    //  fprintf(snp_file,"%i  %i  %c %c  \n",j,i-1,bin2NT[first_unindexed & 3],bin2NT[nt]);
                                    check = true;
                                    break;
                                }
                            }
                            //check other kmers ?  check =false if problem
                            //for higher confidence it would also be possible to test other overlapping kmers

                            if(check)
                            {
                                readseq[ii-1]=bin2NT[nt]; //correc sequence in ram
                                // printf("error found pos %i  read %i\n",ii, _bloocoo._seq_num);
                                _local_nb_errors_corrected ++; //not thread safe
                            }

                        }
                        else if ((tai_not_indexed < sizeKmer) && first_gap && tai_not_indexed>0)
                            // an error at the beginning of the read : tai_not_indexed is smaller even for a single snp
                            //and correct it from the right
                        {
                            check = false;

                            kmer_type kmer_end_rev, tempkmer;
                            kmer_end_rev =revcomp(kmer_end, sizeKmer);
                            //kmer_end.printASCII(sizeKmer);
                            //kmer_end_rev.printASCII(sizeKmer);

                            int strand =1;
                            int nt2;
                            for(nt2=0; nt2<4; nt2++)
                            {
                                tempkmer = model.codeSeedRight (kmer_end_rev, nt2, Data::INTEGER); //to go left : go right with reverse of kmer_end
                                if(_bloom.contains(tempkmer)) //kmer is indexed
                                {
                                    //tempkmer.printASCII(sizeKmer);
                                    check = true;

                                    break;
                                }
                            }


                            //correc error
                            if(check)
                            {
                                //printf("error found pos %i  read %i\n",tai_not_indexed-1, _bloocoo._seq_num);
                                //printf("%c \n",readseq[tai_not_indexed-1]);

                                readseq[tai_not_indexed-1]=bin2NTrev[nt2]; //correc in ram
                                _local_nb_errors_corrected ++; //not thread safe
                                //printf("%c \n",readseq[tai_not_indexed-1]);

                            }
                        }

                    }


                    first_gap = false;

                    //   if(tai_indexed==1) tai_previous_break = tai_not_indexed; // begin of indexed zone    ( if tai_not_indexed ?)
                    if (tai_indexed > 1) // do not reset  tai_not_indexed if a single positive (might be a FP)
                        tai_not_indexed = 0; //reset previous hole size
                    if (tai_indexed==1)// begin of indexed zone
                    {
                        //printf("not indexed %lli\n",tai_not_indexed);
                        tai_not_indexed++; //  treat a single solidkmer  as an erroneous kmer
                    }
                }
                else //kmer contains an error
                {
                    if(tai_indexed==1) //begin of not indexed zone,  previous kmer was an isolated positive, probably a FP
                    {

                    }
                    else if(tai_indexed > 1) // begin of not indexed zone
                    {
                        kmer_begin = previous_kmer ; //kmer_begin is the last kmer indexed before the hole
                        first_unindexed = graine;
                        //printf("indexed %lli\n",tai_indexed);
                    }
                    tai_not_indexed ++;
                    tai_indexed =0;


                    if(ii == (readlen-sizeKmer)) //end of the read, we should treat this  gap here
                        //correc snp with trad method
                    {
                        int nt;
                        check = false;

                        kmer_type temp_kmer;

                        for(nt=0; nt<4; nt++)
                        {
                            temp_kmer = model.codeSeedRight (kmer_begin, nt, Data::INTEGER);
                            if(_bloom.contains(temp_kmer)) //kmer is indexed
                            {
                                check = true;
                                break;
                            }
                        }
                        //verif other kmers :
                        //todo

                        if(check)
                        {

                            //correc error
                            //printf("%c \n",readseq[readlen-1 - tai_not_indexed +1]);

                            readseq[readlen - tai_not_indexed ]=bin2NT[nt]; //correc sequence in ram
                            //printf("%c \n",readseq[readlen-1 - tai_not_indexed +1]);

                            //printf("error found pos %i  read %i\n",readlen-1 - tai_not_indexed +1, _bloocoo._seq_num);
                            _local_nb_errors_corrected ++; //not thread safe
                        }
                    }
                }

                previous_kmer = graine; // should be the kmer in its original sense
            } // end of kmers iteration over the read

            //                if(tai_indexed)
            //                    printf("indexed %lli\n",tai_indexed);
            //                else
            //                    printf("not indexed %lli\n",tai_not_indexed);

        }
#endif

        // TO BE IMPROVED... SHOULD USE A SYNCHRONIZED BANK INSTEAD OF DEALING WITH SYNCHRONIZATION HERE...
        if (_synchro)  { _synchro->lock ();   }
        _outbank.insert (s); //output corrected sequence
        if (_synchro)  { _synchro->unlock (); }

        // printf("%s\n",s.getDataBuffer());
        // printf("%s\n",readseq);

    }

    CorrectReads (Bloom<kmer_type>& bloom, Bank& outbank, Bloocoo & bloocoo, ISynchronizer* synchro, u_int64_t& nb_errors_corrected)
        : _bloom(bloom), _outbank(outbank), _bloocoo(bloocoo), _synchro(synchro),
          _total_nb_errors_corrected (nb_errors_corrected), _local_nb_errors_corrected(0),
          model(_bloocoo._kmerSize), itKmer(model)
    {}

    ~CorrectReads ()
    {
        /** We increase the global number of corrected errors. */
        __sync_fetch_and_add (&_total_nb_errors_corrected, _local_nb_errors_corrected);
    }

    Bloom<kmer_type>& _bloom; // the bloom containing the solid kmers
    Bank&             _outbank; // the bloom containing the solid kmers
    Bloocoo &         _bloocoo;
    ISynchronizer*    _synchro;

    KmerModel           model;
    KmerModel::Iterator itKmer;


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
    _parser->add (new OptionOneParam (DSK::STR_KMER_SIZE,   "size of a kmer",   true));
    _parser->add (new OptionOneParam (DSK::STR_DATABASE,    "database",         true));
    _parser->add (new OptionOneParam (DSK::STR_SOLID_KMERS, "solid kmers file", false));
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
    _kmerSize  = _input->getInt (DSK::STR_KMER_SIZE);
    _solidFile = _input->getStr (DSK::STR_SOLID_KMERS);

    /*************************************************/
    /** We create a bloom with inserted solid kmers. */
    /*************************************************/
    Bloom<kmer_type>* bloom = createBloom ();
    LOCAL (bloom);

    //iterate over initial file
    Bank inbank (_input->getStr(DSK::STR_DATABASE));

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
    string fileName = _input->getStr(DSK::STR_DATABASE) + "_corrected";
    Bank outbank (fileName);

    u_int64_t total_nb_errors_corrected = 0;

    /*************************************************/
    // We iterate over sequences and correct them
    /*************************************************/
    {
        TIME_INFO (_timeInfo, "sequences correction");

        /** We create a shared synchronizer for writing in output file. */
        ISynchronizer* synchro = System::thread().newSynchronizer();

        _dispatcher->iterate (*itSeq,  CorrectReads (*bloom, outbank, *this, synchro, total_nb_errors_corrected));

        delete synchro;
    }

    /*************************************************/
    // We gather some statistics.
    /*************************************************/
    _info->add (1, "result");
    _info->add (2, "nb errors corrected", "%ld", total_nb_errors_corrected);
    _info->add (2, "corrected file",      fileName);
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
    TIME_INFO (_timeInfo, "fill bloom filter");

    double lg2 = log(2);
    float NBITS_PER_KMER = log (16*_kmerSize*(lg2*lg2))/(lg2*lg2);

    u_int64_t solidFileSize = (System::file().getSize(_solidFile) / sizeof (kmer_type));

    u_int64_t estimatedBloomSize = solidFileSize * NBITS_PER_KMER;
    if (estimatedBloomSize ==0 ) { estimatedBloomSize = 1000; }

    /** We create the kmers iterator from the solid file. */
    Iterator<kmer_type>* itKmers = createIterator<kmer_type> (
        new IteratorFile<kmer_type> (_solidFile),
        solidFileSize,
        "fill bloom filter"
    );
    LOCAL (itKmers);

    /** We instantiate the bloom object. */
    BloomBuilder builder (itKmers, estimatedBloomSize, (int)floorf (0.7*NBITS_PER_KMER));
    Bloom<kmer_type>* bloom = builder.build ();

    /** We return the created bloom filter. */
    return bloom;
}

