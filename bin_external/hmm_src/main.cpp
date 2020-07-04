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

float compute_logodds_and_fill_in_seqmap(Pset &pset, Data &data, SequenceMap& seqmap, Settings &settings, SequenceMap::GMS2_GROUP group, int bac_arc) {
    
    seqmap.Init(data.flag, pset.min_gene_length);
    seqmap.CalcGC(data);

//    if (evidence.data.size())
//        seqmap.AddCodingEvidence(data.evi_dir_orf, data.evi_rev_orf);
    
    if (pset.native[0])
        seqmap.CalcStarts(pset.native[0], data.nt, group);
    else
    {
        if (bac_arc == 0) {
            seqmap.CalcStartsGC(pset.first, data.nt, group);
        }
        else
            seqmap.CalcStartsGC(pset.second, data.nt, group);
    }

    if (pset.native[0])
    {
        seqmap.CalcLogP(pset.native, data.nt, NATIVE_TYPE);
        seqmap.CalcLogodd(pset.native, data.nt, NATIVE_TYPE, group);
    }

    if (pset.first[0])
    {
        seqmap.CalcLogP(pset.first, data.nt, ATYPICAL_TYPE_1);
        seqmap.CalcLogodd(pset.first, data.nt, ATYPICAL_TYPE_1, group);
    }

    if (pset.second[0])
    {
        seqmap.CalcLogP(pset.second, data.nt, ATYPICAL_TYPE_2);
        seqmap.CalcLogodd(pset.second, data.nt, ATYPICAL_TYPE_2, group);
    }

    seqmap.Run(settings.hmm.best_start_before_dp, settings.hmm.delta);

    if (pset.native[0])
        seqmap.AddSiteInfo(pset.native[0], data.nt);
    
    if (bac_arc == 0)
        seqmap.AddSiteInfoAtypical(pset.first, data.nt, group);
    else
        seqmap.AddSiteInfoAtypical(pset.second, data.nt, group);
    
    return seqmap.final_logodd;
}

// ----------------------------------------------------
int main( int argc, char** argv )
{
	try
	{
		std::string VERSION = "1.23";

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

				

                // test multiple groups
                SequenceMap::GMS2_GROUP best_group = SequenceMap::NONE;
                char best_label = 'N';
                float best_score = -10000000000;
                SequenceMap::GMS2_GROUP all_groups [] = {
                    SequenceMap::NONE, SequenceMap::A, SequenceMap::B, SequenceMap::C, SequenceMap::D, SequenceMap::X
                };
                char group_labels []  = {'N', 'A', 'B', 'C', 'D', 'X'};
                int best_type = 0;
                
                for (int bac_arc = 0; bac_arc < 2; bac_arc+=1) {
                    for (int group_idx = 0; group_idx < 6; group_idx++) {
                        
                        // bacteria cannot be group D
                        if (bac_arc == 0 && (all_groups[group_idx] == SequenceMap::D))
                            continue;
                        // archaea only groups A and D
                        else if (bac_arc == 1 && (all_groups[group_idx] != SequenceMap::A && all_groups[group_idx] != SequenceMap::D))
                            continue;
                        
                        float current_score = compute_logodds_and_fill_in_seqmap(pset, data, seqmap, settings, all_groups[group_idx], bac_arc);
                        if (current_score > best_score) {
                            best_score = current_score;
                            best_group = all_groups[group_idx];
                            best_label = group_labels[group_idx];
                            best_type = bac_arc;
                        }
                        std::cout << (bac_arc == 0 ? "Bacteria" : "Archaea") << "\t" << group_labels[group_idx] << "\t" << current_score << std::endl;
                    }
                }
                
                std::cout << (best_type == 0 ? "Bacteria" : "Archaea") << "\t" << "Best group: " << best_label << "\t" << best_score << std::endl;
                
                // rerun with best group
                compute_logodds_and_fill_in_seqmap(pset, data, seqmap, settings, best_group, best_type);
				

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

