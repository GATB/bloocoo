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

#ifndef __bloocoo_ReadWriter__
#define __bloocoo_ReadWriter__

/********************************************************************************/

#include <gatb/tools/misc/impl/Tool.hpp>
#include <gatb/tools/math/Integer.hpp>
#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>
#include <gatb/tools/designpattern/impl/Command.hpp>
#include <gatb/tools/misc/impl/TimeInfo.hpp>
#include <gatb/tools/misc/impl/OptionsParser.hpp>
#include <gatb/bank/impl/Banks.hpp>
#include <gatb/kmer/impl/Model.hpp>
#include <gatb/tools/collections/impl/Bloom.hpp>
#include <gatb/tools/collections/impl/BagFile.hpp>
#include <gatb/tools/collections/impl/BagCache.hpp>
#include <string>
#include <sstream>
#include <iostream>
#include <pthread.h>


/********************************************************************************/
/** NOTE: we should not include namespaces here => only to make user life easier... */
using namespace gatb::core;
using namespace gatb::core::tools;
using namespace gatb::core::bank;
using namespace gatb::core::kmer::impl;
using namespace gatb::core::system;
using namespace gatb::core::system::impl;
using namespace gatb::core::bank::impl;
using namespace gatb::core::tools::collections;

/********************************************************************************/

class OrderedBankWriter;

typedef struct
{
    OrderedBankWriter * obw;
} thread_args;

void * writer(void * args);





class OrderedBankWriter : public SmartPointer
{
public:
    
    /** Constructor. */
    OrderedBankWriter (Bag<Sequence> * bank, size_t buffsize)
    : _writer_available(1), _buffer_full(0),_to_be_written(0), _buffWrite (buffsize), _buffReceive (buffsize),_ref(0),  _nbMax(buffsize), _idx(0), _base(0),_writer_kill_order(0)
    {
        setRef(bank);
        
        //create writer thread
        t_arg.obw = this;
        
       // printf ("  t_arg.obw %p    bank %p ",t_arg.obw,_ref);
        pthread_mutex_init(&writer_mutex, NULL);
        pthread_cond_init (&writer_available_cond, NULL);

        pthread_cond_init (&buffer_full_cond, NULL);
        
        pthread_create (&_thread, NULL,  writer, &t_arg);
	 // printf("--OrderedBankWriter constructor, create thread  %p  ,    this %p    cond %p   targs %p--\n",&_thread,this,&buffer_full_cond, &t_arg );

    }
    
    /** Destructor. */
    ~OrderedBankWriter ()
	{
		//this should not be called before  writer thread has been killed !!
		
		
		//////////// first wait for writer thread to finish
		pthread_join(_thread,NULL);
		
		//ici mettre code attente threads
		//printf(" OrderedBankWriter destructor,   this %p   \n",this );
		
		pthread_mutex_destroy(&writer_mutex);

        setRef(0);
    }
    
    
    void insert(const Sequence&  seq ) ;
    
    void incDone (int nbdone);

    void waitForWriter ();
    void FlushWriter ();

    
    thread_args t_arg;
    pthread_mutex_t writer_mutex;
    pthread_cond_t writer_available_cond;

    pthread_cond_t buffer_full_cond;
    int _writer_available ;
	int _writer_kill_order ;

    int _buffer_full;
    int _to_be_written;
    
    std::vector<Sequence >      _buffWrite;
    std::vector<Sequence >      _buffReceive;


    

    
    Bag<Sequence> *        _ref;

protected:
    
    
    void setRef (Bag<Sequence>* ref)
    {
        _ref=ref;
    }



    //they will be swapped _buffWrite.swap(_buffReceive);
    



    size_t      _nbMax; // buffer size
    size_t      _idx;   //number of sequences stored in the buffer
    size_t      _base;
    
    pthread_t  _thread; // the writer thread


    
private:
    //assign
    OrderedBankWriter& operator=(const OrderedBankWriter& bk){
        return *this;
    }
    //copy construct
    OrderedBankWriter(const OrderedBankWriter& bk){ }
	
    
};




#endif /* defined(__bloocoo_xcode__ReadWriter__) */
