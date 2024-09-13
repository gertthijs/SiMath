/*
 *  Notifier.h
 *
 *  Created by Gert Thijs on 11/04/06.
 *  Copyright 2006 Silicos NV. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *     + Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     + Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY SILICOS NV AND CONTRIBUTORS ``AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL SILICOS NV AND CONTRIBUTORS BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
		* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
 *
 */

#ifndef _SIMATH_NOTIFIER_H__
#define _SIMATH_NOTIFIER_H__

#include <iostream>

namespace SiMath
{
	//------------------------------------------------------------------------------
	/**
	 * \class Notifier
	 * \brief Class providing a general framework of notifying progress and log
	 *        messages.
	 *
	 * This simple class provides an interface for notifying the outside world of
	 * a class that something was changed in that class. Every class (or function)
	 * can support this type of communication by using this general notifier. 
	 * The user of that class (or function) can later specify its own notifiers 
	 * without changing the implementation of that class.
	 *
	 * \code
	 *   void aFunction(Notifier* n=0) {
	 *     for(int i=0 ; i<100 ; ++i) {
	 *       // do something
	 *       if(n!=0)
	 *       {
	 *         n->log("iteration completed without problem.");
	 *         n->progress(i);
	 *       }
	 *     }  
	 *   }
   *
   *   class EmailNotifier: public Notifier { 
	 *     public:
	 *       void progress(int i) {
	 *         if(i%10 == 0)
	 *           sendEmail("Notification!", ...);
	 *       }
	 *       void log(string s) {
	 *         // do nothing
	 *       }
	 *   }
   *
   *   EmailNotifier* n = new EmailNotifier;
   *   aFunction(n);
   * \endcode
   *
	 * As we take a look to the previous code example, we can see that the
	 * function \c aFunction caused an email notification every tenth generation
	 * without even knowing what an email is. It is just able to notify 
	 * 'something'.
	 *
	 * This approach can be very useful if we want to incorporate code into a GUI.
	 * In the previous code example we simply make a new notifier that corresponds
	 * with some kind of Progress Bar and our original code still remains 
	 * independent of any GUI at all.
	 *
	 * If a class (or function) provides notifier(s), it is recommended to mention
	 * it in the documenation, so that it will be easier for further users to 
	 * define their own notifiers.
	 *
	 * \author Jonatan Taminau
	 */
	class Notifier
	{
		public:
		
			/**
			 * \brief Notification of progress. Default written to \c std::cerr.
			 */
			virtual void progress(int i) {
				std::cerr << i << std::endl;
			}
		
			/**
			 * \brief Notification of log messages. Default written to \c std::cerr.
			 */
			virtual void log(std::string s) {
				std::cerr << s << std:endl;
			}
	};

};

#endif _SIMATH_NOTIFIER_H__
