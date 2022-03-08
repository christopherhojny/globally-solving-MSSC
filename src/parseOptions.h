// -*- C++ -*-
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*    This file is part of the program colorbitopt                           */
/*                                                                           */
/*    an implementation of a branch-and-cut algorithm to solve the           */
/*    coloring problem by symmetry breaking methods based on orbitopes.      */
/*                                                                           */
/*    Copyright (C) 2005-2014  Marc Pfetsch                                  */
/*                                                                           */
/*                                                                           */
/*    Based on SCIP  --- Solving Constraint Integer Programs                 */
/*                                                                           */
/*    Copyright (C) 2002-2014 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*       mailto: scip@zib.de                                                 */
/*       SCIP is distributed under the terms of the SCIP Academic Licence,   */
/*       see file COPYING in the SCIP distribution.                          */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   parseOptions.h
 * @brief  Function to obtain options from a command line interface string
 * @author Marc Pfetsch
 *
 * There are many tools/libraries that allow the parsing of options as they are given on a command
 * line:
 *
 * - getopt (C)
 * - argvp (C)
 * - parse_options (C++_
 *
 * These libraries are, however, too complex for many purposes. One simply wants to say how many
 * required parameters there are and what are the optional parameters. This small tool should fill
 * the gap.
 *
 * This tools steals code and concepts from the corresponding "argv_option" in polymake.
 */

#ifndef PARSEOPTIONS_H
#define PARSEOPTIONS_H

#include <vector>
#include <cstring>


/**@struct parseOption
 * @brief Store information for options
 *
 * This struct stores an option like '-output <...>'. @c takesValues says whether the option needs an argument or not.
 */
struct parseOption
{
   const char*           name;               //!< name of option
   bool                  takesValue;         //!< whether the option requires a value to be given
   const char*           value;              //!< string possibly containing the value

   operator const char* () const
   {
      return value;
   }
};



/** Parses the options given in @c argv and stores them in @c options
 *  @returns '@c true' if succesful
 *
 *  Options start with '-'.
 *
 *  @note We do not change argv or argc.
 */
bool extractOptions(
   int                   argc,               //!< number of arguments
   const char*const*     argv,               //!< array of arguments
   int                   nOptions,           //!< number of options (fixed in advance)
   parseOption*          options,            //!< array of options (fixed in advance)
   int                   minOtherArgs,       //!< mininum number of other arguments
   int                   maxOtherArgs,       //!< maximum number of other arguments
   std::vector<std::string>& otherArgs       //!< array of other arguments (on output)
   )
{  /*lint --e{850}*/
   otherArgs.clear();
   if ( argc < minOtherArgs )
      return false;

   int nOtherArgs = 0;
   for (int i = 1; i < argc; ++i)
   {
      // if the argument is not an option
      if ( argv[i][0] != '-' )
      {
	 otherArgs.push_back(std::string(argv[i]));
	 ++nOtherArgs;
	 if ( nOtherArgs > maxOtherArgs )
	    return false;
	 continue;
      }
      // check whether we recognize the option
      bool recognized = false;
      for (int j = 0; j < nOptions; ++j)
      {
	 if (! std::strcmp(argv[i], options[j].name))
	 {
	    if (options[j].takesValue)
	    {
	       if (i == argc-1)
		  return false;
	       options[j].value = argv[++i];
	    }
	    else
	       options[j].value = "yes";
	    recognized = true;
	    break;
	 }
      }
      if (!recognized)
	 return false;
   }
   if ( nOtherArgs < minOtherArgs )
      return false;
   return true;
}



/** Parses the options given in @c argv and stores them in @c options
 *  @returns '@c true' if succesful
 *
 *  For this variant, we do not need to specify the number of options.
 */
template <int nOptions> inline
bool extractOptions(
   int                   argc,               //!< number of arguments
   const char**          argv,               //!< array of arguments
   parseOption(&options)[nOptions],          //!< array of options (fixed in advance)
   int                   minOtherArgs,       //!< mininum number of other arguments
   int                   maxOtherArgs,       //!< maximum number of other arguments
   std::vector<std::string>& otherArgs       //!< array of other arguments (on output)
   )
{
   return extractOptions(argc, argv, nOptions, options, minOtherArgs, maxOtherArgs, otherArgs);
}

#endif
