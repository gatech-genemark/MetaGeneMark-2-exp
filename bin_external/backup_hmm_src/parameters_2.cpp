// ====================================================
// Author: Alex Lomsadze
// Copyright: Georgia Institute of Technology
// Copyright: GeneProbe Inc.
// File: parameters_2.cpp
// Project: GeneMark.hmm-2 (no introns)
// ====================================================

#include <string>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <utility>

using std::string;
using std::map;
using std::cout;
using std::endl;
using std::ifstream;
using std::stringstream;
using std::pair;

#include "parameters_2.h"
#include "parameter_parser_2.h"
#include "exit.h"

// ----------------------------------------------------
Parameters::Parameters( Settings const & settings, Logger * const logger ) : logger(logger)
{
	if ( logger->verbose ) logger->Print( 0, "# Loading parameters ..." );

	native_filename = settings.hmm.native_filename;
	mgm_filename    = settings.hmm.mgm_filename;
	tis_filename    = settings.hmm.tis_filename;

	// the same parameter can be specified in three locations
	// location processing order: defaults < from model file < from command line

	to_native = 0.5;
	to_mgm = 0.5;
	to_atypical_first  = 0.5;
	to_atypical_second = 0.5;

	// load from file
	Load();

	// load from command line
	if ( ! settings.hmm.bac_prob_defaulted )  to_atypical_first  = settings.hmm.bac_prob;
	if ( ! settings.hmm.arc_prob_defaulted )  to_atypical_second = settings.hmm.arc_prob;
	if ( ! settings.hmm.nat_prob_defaulted )  to_native          = settings.hmm.native_prob;
	if ( ! settings.hmm.mgm_prob_defaulted )  to_mgm             = settings.hmm.mgm_prob;

	if ( logger->debug )   logger->Print( 10, Summary() );
	if ( logger->verbose ) logger->Print( 0, "# Loading parameters ... done" );
}
// ----------------------------------------------------
void Parameters::Load(void)
{
	ParameterParser  parser;

	if ( !mgm_filename.empty() )
	{
		parser.PutFileToBuffer( mgm_filename );

		// "__B_GC" and "__A_GC" are expected

		LoadSection( parser, parser.buffer, "__" );
	}

	if ( !native_filename.empty() )
	{
		parser.PutFileToBuffer( native_filename );

		LoadSection( parser, parser.buffer, "__NATIVE" );
	}

	if ( !tis_filename.empty() )
	{
		parser.PutFileToBuffer( tis_filename );

		LoadSection( parser, parser.buffer, "__START" );
	}

	if ( !model_set.size() )
		Exit( "error, model set is empty" );
}
// ----------------------------------------------------
void Parameters::LoadSection( ParameterParser & parser, std::string & buffer,  const std::string label )
{
	// for logger
	unsigned int set_size = model_set.size();

	// Parse 'buffer' into sections with "__" separator

	parameter_map  section;

	parser.LoadFromString( section, buffer, "__" );

	parameter_map::iterator itr     = section.begin();
	parameter_map::iterator itr_end = section.end();

	for( ; itr != itr_end; ++itr )
	{
		string section_name = itr->first;

		// skip sections whithout partial match to the label

		if ( section_name.compare( 0, label.size(), label ) ) 
			continue;

		if (logger->debug) logger->Print( 10, "Loading section: ", section_name  );

		// put empty model into map
	
		if ( model_set.find( section_name ) == model_set.end() )
		{
			model_set.insert( pair<string, Model>( section_name, Model(logger) ) );
		}
		else
			Exit( "error, section duplication found in parameter file:", section_name );

		// Parse model into variables using '$' separator

		parameter_map  par;

		parser.LoadFromSubString( par, buffer, "$", itr->second.first, itr->second.second );

		map< string, Model >::iterator current = model_set.find( section_name );

		// Fill in empty model

		// =====================================================================================================
		{
			if (parser.IsKey(par, "$NAME"))            current->second.model_name         = parser.asString ( par, "$NAME", buffer );
			if (parser.IsKey(par, "$GCODE"))           current->second.gcode              = parser.asPInt   ( par, "$GCODE", buffer );
			if (parser.IsKey(par, "$GENE_MIN_LENGTH")) current->second.gene_min_length    = parser.asPInt   ( par, "$GENE_MIN_LENGTH", buffer );
			if (parser.IsKey(par, "$COD_ORDER"))       current->second.order_cod          = parser.asPInt   ( par, "$COD_ORDER", buffer );
			if (parser.IsKey(par, "$NON_ORDER"))       current->second.order_non          = parser.asPInt   ( par, "$NON_ORDER", buffer );
			if (parser.IsKey(par, "$COD_P_N")) current->second.probability_N_in_coding    = parser.asPDouble( par, "$COD_P_N", buffer );
			if (parser.IsKey(par, "$NON_P_N")) current->second.probability_N_in_noncoding = parser.asPDouble( par, "$NON_P_N", buffer );
			if (parser.IsKey(par, "$NON_DURATION_DECAY")) current->second.noncoding_duration_decay = parser.asPDouble( par, "$NON_DURATION_DECAY", buffer );
			if (parser.IsKey(par, "$COD_DURATION_DECAY")) current->second.coding_duration_decay    = parser.asPDouble( par, "$COD_DURATION_DECAY", buffer );

			if (parser.IsKey(par, "$BUILD"))   current->second.build = parser.asString(par, "$BUILD", buffer);

			if ( parser.IsKey( par, "$ATG" ) )  current->second.pATG = parser.asPDouble( par, "$ATG", buffer );
			if ( parser.IsKey( par, "$GTG" ) )  current->second.pGTG = parser.asPDouble( par, "$GTG", buffer );
			if ( parser.IsKey( par, "$TTG" ) )  current->second.pTTG = parser.asPDouble( par, "$TTG", buffer );
			if ( parser.IsKey( par, "$TAA" ) )  current->second.pTAA = parser.asPDouble( par, "$TAA", buffer );
			if ( parser.IsKey( par, "$TAG" ) )  current->second.pTAG = parser.asPDouble( par, "$TAG", buffer );
			if ( parser.IsKey( par, "$TGA" ) )  current->second.pTGA = parser.asPDouble( par, "$TGA", buffer );
		
			if (parser.IsKey(par, "$NON_MAT") || parser.IsKey(par, "$COD_MAT"))
			{
				current->second.ReserveSpace();

				parser.asVectorOfDoublesWithLabel(par, "$NON_MAT", buffer, current->second.non);
				parser.as3VectorOfDoublesWithLabel(par, "$COD_MAT", buffer, current->second.cod1, current->second.cod2, current->second.cod3);
			}
		}
		// =====================================================================================================
		{
			LoadSite( parser, par, buffer, "$RBS",          &current->second.RBS,                  true );
			LoadSite( parser, par, buffer, "$PROMOTER",     &current->second.Promoter,             true );
			LoadSite( parser, par, buffer, "$SC",           &current->second.StartContent,         false );
			LoadSite( parser, par, buffer, "$SC_RBS",       &current->second.StartContentRBS,      false );
			LoadSite( parser, par, buffer, "$SC_PROMOTER",  &current->second.StartContentPromoter, false );
		}
		// =====================================================================================================
		{
			if ( parser.IsKey( par, "$TO_ATYPICAL_FIRST_BACTERIA" ) )  to_atypical_first  = parser.asPDouble( par, "$TO_ATYPICAL_FIRST_BACTERIA", buffer );
			if ( parser.IsKey( par, "$TO_ATYPICAL_SECOND_ARCHAEA" ) )  to_atypical_second = parser.asPDouble( par, "$TO_ATYPICAL_SECOND_ARCHAEA", buffer );
			if ( parser.IsKey( par, "$TO_NATIVE" ) )                   to_native          = parser.asPDouble( par, "$TO_NATIVE", buffer );
			if ( parser.IsKey( par, "$TO_MGM" ) )                      to_mgm             = parser.asPDouble( par, "$TO_MGM", buffer );

			if ( parser.IsKey( par, "$GENE_GC_DIST" ) && parser.asBool(par, "$GENE_GC_DIST", buffer ) )
			{
				to_native_by_gc.assign( 101, 0 );

				parser.asVectorOfDoublesWithPos( par, "$GENE_GC_DIST_MAT", buffer, to_native_by_gc );
			}
		}
	}

	if (logger->debug) logger->Print( 10, "Sections loaded:", label,  model_set.size() - set_size  );
}
// ----------------------------------------------------
void Parameters::LoadSite( ParameterParser & parser, parameter_map & par, std::string & buffer, std::string const label, Site * ptr, bool with_dur )
{
	if ( parser.IsKey( par, label ) && parser.asBool( par, label, buffer ) )
	{
		ptr->order    = parser.asPInt( par, label + "_ORDER",   buffer );
		ptr->width    = parser.asPInt( par, label + "_WIDTH",   buffer );
		ptr->margin   = parser.asInt ( par, label + "_MARGIN",  buffer );

		if ( with_dur )
			ptr->max_dur  = parser.asPInt( par, label + "_MAX_DUR", buffer ) + 1;

		ptr->ReserveSpace();

		parser.asMatrixOfDoublesWithLabel( par, label + "_MAT", buffer, ptr->matrix );

		if ( with_dur)
			parser.asVectorOfDoublesWithPos( par, label + "_POS_DISTR", buffer, ptr->duration );

		ptr->is_valid = true;
	}
}
// ----------------------------------------------------
string Parameters::Summary(void)
{
	stringstream s;

	s << "# START Parameters:" << endl;
	s << "parameter sets loaded into 'model_set' : " << model_set.size() << endl;

	s << "to_native: " << to_native << endl;
	s << "to_mgm: " << to_mgm << endl;
	s << "to_atypical_first: " << to_atypical_first << endl;
	s << "to_atypical_second: " << to_atypical_second << endl;

	s << "# STAR KL of model" << endl;

	map< string, Model >::iterator  itr     = model_set.begin();
	map< string, Model >::iterator  itr_end = model_set.end();

	for( ; itr != itr_end; ++itr )
	{
		s << itr->first << " " << itr->second.KL() << endl;
	}

	s << "# END KL of model" << endl;
	s << "# END Parameters " << endl;

	return s.str();
}
// ----------------------------------------------------
void Parameters::Initialize( std::map< std::string, Model > & target )
{
	map< string, Model >::iterator  itr     = target.begin();
	map< string, Model >::iterator  itr_end = target.end();

	for( ; itr != itr_end; ++itr )
	{
		itr->second.Initialize();
	}
}
// ----------------------------------------------------
