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


#include "ReadWriter.h"


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

//        __sync_fetch_and_add (&_total_nb_errors_corrected, _local_nb_errors_corrected);

void OrderedBankWriter::insert(const Sequence&  seq)
{
    //seq id begins at 1
    // or 0 ?
    size_t id;  
   // Sequence * tempseq;

    while (1)
    {
         id  = seq.getIndex() - _base;
        if(id < _nbMax )
        {
            
#ifdef ASSIGN
            _buffReceive[id] =  seq;

#else
            _buffReceive[id] =  (new Sequence(seq));

#endif
            
           // printf("inserting into %i    ::  %p\n",id,_buffReceive[id]);
         //   tempseq= new Sequence(seq);
           // _buffReceive[id] =  *tempseq;
          //  delete tempseq;
            break;
        }
        else
        {
            //the current thread has finished the block and is now ahead of others,
            //wait for other threads to finish filling the buffer

            
            pthread_mutex_lock(&writer_mutex);
           // while (_buffer_full==0)
            while (id >= _nbMax )
            {
                pthread_cond_wait(&buffer_full_cond, &writer_mutex);
                 id  = seq.getIndex() - _base;
            }
            pthread_mutex_unlock(&writer_mutex);
            
    
        }
    }
    
    
}


void OrderedBankWriter::incDone (int nbdone) 
{
    
    size_t old_val;


    old_val = __sync_fetch_and_add (& _idx, nbdone);
   // printf("\n ++ ninc oldval  %i base = %i ++ \n",old_val,_base);

    // buffer completely filled : activates writer thread
    if( (old_val + nbdone) ==_nbMax)  
    {
       // printf("\n ++ ninc done %zd  and buff done, base = %zd ++ \n",old_val,_base);

        ///// wait for writer to be free
        pthread_mutex_lock(&writer_mutex);
        while (_writer_available==0) {
            pthread_cond_wait(&writer_available_cond, &writer_mutex);
        }
        //buffReceive is full, will be written

        _buffWrite.swap(_buffReceive);
        
        _to_be_written = _nbMax;
        _idx=0;
        _base += _nbMax;
        
        //signal writer he should write buffer
        _buffer_full = 1;
        pthread_cond_broadcast(&buffer_full_cond);
        pthread_mutex_unlock(&writer_mutex);


        
    }
    
}

void OrderedBankWriter::waitForWriter ()
{
    pthread_mutex_lock(&writer_mutex);
    while (_writer_available==0) {
   //     printf("worker going to sleep, waiting for writer to be finished .. \n");
        pthread_cond_wait(&writer_available_cond, &writer_mutex);
    }
    pthread_mutex_unlock(&writer_mutex);
 //   printf("end waitForWriter \n");

}


//guaranteed by design only one thread at the end will call this
void OrderedBankWriter::FlushWriter ()
{

    if( _idx)
    {
      //  printf("\n flushing %zd base = %zd ++ \n",_idx,_base);
        
        ///// wait for writer to be free
        pthread_mutex_lock(&writer_mutex);
        while (_writer_available==0) {
            pthread_cond_wait(&writer_available_cond, &writer_mutex);
        }
        //buffReceive is full, will be written
        
        _buffWrite.swap(_buffReceive);
        _to_be_written = _idx;
        
        //signal writer he should write buffer
        _buffer_full = 1;
        pthread_cond_broadcast(&buffer_full_cond);
        pthread_mutex_unlock(&writer_mutex);
    
    }
    
    pthread_mutex_lock(&writer_mutex);
    while (_buffer_full==1) {
        pthread_cond_wait(&writer_available_cond, &writer_mutex);
    }
    pthread_mutex_unlock(&writer_mutex);
}

// writer thread gets a pointer to  orderbankwriter object
void * writer(void * args)
{

    thread_args *targ = (thread_args*) args;
    OrderedBankWriter * obw = targ->obw;
    Bank * outbank = obw->_ref;
    int  * to_be_written =  &(obw->_to_be_written);
#ifdef ASSIGN
    std::vector<Sequence> * _buffWrite = &(obw->_buffWrite);
#else
    std::vector<Sequence*> * _buffWrite = &(obw->_buffWrite);

#endif
    


    

    //when waken, writes content of  vector  _buffWrite to bank
    while(1)
    {
        pthread_mutex_lock(&(obw->writer_mutex));
        while (obw->_buffer_full==0)
        {
           // printf(" writer thread  going to sleep .. obw %p   cond %p\n",obw,&(obw->buffer_full_cond));
            pthread_cond_wait(&(obw->buffer_full_cond), &(obw->writer_mutex));
        }
        obw->_writer_available = 0;
      //  printf(" writer thread  awaken !! ..  will write %i elems \n",*to_be_written);

        pthread_mutex_unlock(&(obw->writer_mutex));
        

        //writes the buffer
#ifdef ASSIGN
        for ( std::vector<Sequence>::iterator it = _buffWrite->begin(); (it != _buffWrite->end()) && (*to_be_written) ; it++)
#else
        for ( std::vector<Sequence*>::iterator it = _buffWrite->begin(); (it != _buffWrite->end()) && (*to_be_written); it++)
#endif
  
        {

#ifdef ASSIGN
             outbank->insert((*it));
#else
                outbank->insert(*(*it));
#endif
            //printf("to be written %i \n",*to_be_written);
           ( *to_be_written) --;
        }
        
        //signal buffer has been written
        pthread_mutex_lock(&(obw->writer_mutex));
        obw->_buffer_full=0;
        obw->_writer_available=1;
      //  printf(" writer thread  setting  _buffer_full to 0 .. %i  tobewrit %i \n",obw->_buffer_full,*to_be_written);
        pthread_cond_signal(&(obw->writer_available_cond));
        pthread_mutex_unlock(&(obw->writer_mutex));
    }
    
    
}

