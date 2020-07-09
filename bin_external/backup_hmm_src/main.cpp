// ====================================================
// Author: Alex Lomsadze
// 
// Copyright: Georgia Institute of Technology
// Copyright: GeneProbe Inc.
// 
// Contact information:
//     Mark Borodovsky  borodovsky@gatech.edu
//     Alex Lomsadze    alexl@gatech.edu
// 
// GeneMark.hmm-2 algorithm version 1.0 was published in
// Genome Research 2018
// https://www.ncbi.nlm.nih.gov/pubmed/29773659
// 
// "Modeling leaderless transcription and atypical genes
// results in more accurate gene prediction in prokaryotes."
// 
// by Alexandre Lomsadze, Karl Gemayel, Shiyuyun Tang and Mark Borodovsky
// 
// Project: GeneMark.hmm-2 (no introns)
// File: main.cpp
//
// Projects VERSION is set by first declaration in the main() function
// ====================================================

#include <string>
#include <iostream>

#include "exit.h"
#include "logger.h"
#include "settings_2.h"
#include "sequence_file_2.h"
#include "parameters_2.h"
#include "pset_2.h"
#include "data_2.h"
#include "sequence_map_2.h"
#include "evidence_2.h"
#include "model_2.h"
#include "output_2.h"

#ifdef LICENSE
#include "check.h"
#endif

// ----------------------------------------------------
int main( int argc, char** argv )
{
	try
	{
		std::string VERSION = "1.21";

#ifdef LICENSE
		VERSION += "_lic";
#endif

		Logger        logger;
		Settings      settings( argc, argv, &logger, VERSION );

#ifdef LICENSE
		// to do:
		// move to new version of OS independent key implementation
		// function should return true/false and support verbose mode
		// function should be called from settings after parsing of parameters and before usage

		char path_to_key[] = "";
		check_timekey(path_to_key);
#endif

		SequenceFile  sequence( settings, &logger );
		Evidence      evidence( settings, &logger );

		if ( sequence.AllAtOnce() )
			evidence.SyncWithSequence( sequence.data );

		Parameters    parameters( settings, &logger );
		Pset          pset( settings, parameters, &logger );

		parameters.Initialize( parameters.model_set );

		Data          data( settings, pset, &logger );
		SequenceMap   seqmap( settings, &logger );
		Output        output( settings, &logger, VERSION );

		output.evi_dir_orf = & data.evi_dir_orf;
		output.evi_rev_orf = & data.evi_rev_orf;

		output.Header1(pset);

		bool run_on_all_contigs = true;
		if (!settings.in.run_on_records.empty() || !settings.in.ignore_records.empty())
			run_on_all_contigs = false;

		// progress bar variables
		int records_in = 0;
		int current_record = 0;
		int progress_bar_at = 0;

		if (settings.progress_bar)
			records_in = sequence.RecordsIn();

		for( FastaVectorItr itr = sequence.data.begin() ; itr != sequence.data.end(); itr = sequence.Next( itr ) )
		{
			if (settings.progress_bar)
			{
				if (current_record == progress_bar_at)
				{
					int progress_percent = int(100 * progress_bar_at / records_in);
					std::cout << "Progress:" << progress_percent << "%" << std::endl;
					progress_bar_at += int(records_in / 10);
				}

				++current_record;
			}

			if (!run_on_all_contigs)
			{
				if (!sequence.RunOnThis(itr->name))
					continue;

				if (sequence.IgnoreThis(itr->name))
					continue;
			}

			output.Header2( itr );

			// skip short sequences
			if (itr->data.size() < pset.min_gene_length)
				continue;

			do
			{
				data.Set(itr->data);

				if (evidence.data.size())
					data.ApplyEvidence(evidence.data, itr->name);

				seqmap.Init(data.flag, pset.min_gene_length);
				seqmap.CalcGC(data);

				if (evidence.data.size())
					seqmap.AddCodingEvidence(data.evi_dir_orf, data.evi_rev_orf);

				if (pset.native[0])
					seqmap.CalcStarts(pset.native[0], data.nt);
				else
				{
					seqmap.CalcStartsGC(pset.first, data.nt);
				}

				if (pset.native[0])
				{
					seqmap.CalcLogP(pset.native, data.nt, NATIVE_TYPE);
					seqmap.CalcLogodd(pset.native, data.nt, NATIVE_TYPE);
				}

				if (pset.first[0])
				{
					seqmap.CalcLogP(pset.first, data.nt, ATYPICAL_TYPE_1);
					seqmap.CalcLogodd(pset.first, data.nt, ATYPICAL_TYPE_1);
				}

				if (pset.second[0])
				{
					seqmap.CalcLogP(pset.second, data.nt, ATYPICAL_TYPE_2);
					seqmap.CalcLogodd(pset.second, data.nt, ATYPICAL_TYPE_2);
				}

				seqmap.Run(settings.hmm.best_start_before_dp, settings.hmm.delta);

				if (pset.native[0])
					seqmap.AddSiteInfo(pset.native[0], data.nt);

			} while (0);

			output.PrintGenes(seqmap.predictions, itr, pset.genetic_code);
			output.Stat(itr, seqmap.final_logodd, seqmap.predictions);
		}

		output.Footer();

		logger.Print( 0, "# GeneMark.hmm ... done" );
	}
	catch( std::bad_alloc const & )
	{
		Exit("error, out of memory");
	}
	catch( std::exception const & e )
	{
		Exit("error,", e.what());
	}

	return 0;
};
// ----------------------------------------------------

