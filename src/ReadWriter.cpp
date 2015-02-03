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


#include <gatb/system/impl/System.hpp>

#include <gatb/bank/impl/Banks.hpp>
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


void OrderedBankWriter::insert(const Sequence&  seq)
{

    size_t id;  

    while (1)
    {
         id  = seq.getIndex() - _base;
        if(id < _nbMax )
        {
        
            _buffReceive[id] =  seq;
            break;
        }
        else
        {
            //the current thread has finished the block and is now ahead of others,
            //wait for other threads to finish filling the buffer
            
            pthread_mutex_lock(&writer_mutex);
			id  = seq.getIndex() - _base;

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
	pthread_mutex_lock(&writer_mutex);

    if( _idx)
    {
      //  printf("\n flushing %zd base = %zd ++ \n",_idx,_base);
        ///// wait for writer to be free
        while (_writer_available==0) {
            pthread_cond_wait(&writer_available_cond, &writer_mutex);
        }

        
        _buffWrite.swap(_buffReceive);
        _to_be_written = _idx;
        
        //signal writer he should write buffer
        _buffer_full = 1;
        pthread_cond_broadcast(&buffer_full_cond);
		
    }
	
	pthread_mutex_unlock(&writer_mutex);

    //wait again for writer to finish flushing
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
    Bag<Sequence>* outbank = obw->_ref;
    int  * to_be_written =  &(obw->_to_be_written);
    std::vector<Sequence> * _buffWrite = &(obw->_buffWrite);

    

    //when waken, writes content of  vector  _buffWrite to bank
    while(1)
    {
        pthread_mutex_lock(&(obw->writer_mutex));
        while (obw->_buffer_full==0)
        {
            pthread_cond_wait(&(obw->buffer_full_cond), &(obw->writer_mutex));
        }
        obw->_writer_available = 0; //writer becomes busy
      //  printf(" writer thread  awaken !! ..  will write %i elems \n",*to_be_written);
        pthread_mutex_unlock(&(obw->writer_mutex));
        

        //writes the buffer
        for ( std::vector<Sequence>::iterator it = _buffWrite->begin(); (it != _buffWrite->end()) && (*to_be_written) ; it++)
        {
            outbank->insert((*it));
            (*to_be_written) --;
        }
        
        
        //signal to others buffer has been written and writer is available again
        pthread_mutex_lock(&(obw->writer_mutex));
        obw->_buffer_full=0;
        obw->_writer_available=1;
      //  printf(" writer thread  setting  _buffer_full to 0 .. %i  tobewrit %i \n",obw->_buffer_full,*to_be_written);
        pthread_cond_signal(&(obw->writer_available_cond));
        pthread_mutex_unlock(&(obw->writer_mutex));
    }
    
    
}

