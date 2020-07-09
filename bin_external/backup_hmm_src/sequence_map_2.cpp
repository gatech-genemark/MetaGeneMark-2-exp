// ====================================================
// Author: Alex Lomsadze
// Copyright: Georgia Institute of Technology
// Release: 2017
// File: sequence_map_2.cpp
// Project: GeneMark.hmm-2 (no introns)
// ====================================================

#include <cassert>
#include <cmath>

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <sstream>

using std::vector;
using std::cerr;
using std::endl;
using std::cout;
using std::stringstream;
using std::string;
using std::map;

#include "sequence_map_2.h"
#include "index_2.h"
#include "data_2.h"
#include "exit.h"

// ----------------------------------------------------
SequenceMap::SequenceMap( Settings & settings, Logger * const logger ) : logger(logger)
{
	overlap_penalty = settings.hmm.overlap_penalty;
	noncoding_state_initiation_and_termination_probability = settings.hmm.noncoding_state_initiation_and_termination_probability;
	final_logodd = LOG_ZERO;
}
// ----------------------------------------------------
void SequenceMap::AddIniTermProb(double p_ini_non)
{
	double logodd_non = LOG_ZERO;
	double logodd_cod = LOG_ZERO;

	if ( p_ini_non == 0 )
	{
		logodd_non = LOG_ZERO;
		logodd_cod = 0;
	}
	else if ( p_ini_non == 1 )
	{
		logodd_non = 0;
		logodd_cod = LOG_ZERO;
	}
	else if ( p_ini_non > 0 && p_ini_non < 1 )
	{
		logodd_non = 0;
		logodd_cod = log( ( 1 - p_ini_non )/p_ini_non );
	}
	else
	{
		cerr << "error, 'SequenceMap::AddIniTermProb()'" << endl;
		exit(1);
	}

	// initiation sites occupy first positions in the map

	for( vector< MapValue >::iterator itr = data.begin(); itr != data.end(); ++itr )
	{
		if (( itr->status & iniCod )||(itr->status & iniNon))
		{
			if ( itr->is_good &&( itr->status & iniCod ))  itr->prob += logodd_cod;
			if ( itr->is_good &&( itr->status & iniNon ))  itr->prob = logodd_non;
		}
		else
			break;
	}

	// termination sites occupy last positions in the map

	for( vector< MapValue >::reverse_iterator itr = data.rbegin(); itr != data.rend(); ++itr )
	{
		if (( itr->status & terCod )||(itr->status & terNon))
		{
			if ( itr->is_good &&( itr->status & terCod ))  itr->prob += logodd_cod;
			if ( itr->is_good &&( itr->status & terNon ))  itr->prob = logodd_non;
		}
		else
			break;
	}
}
// ----------------------------------------------------
void SequenceMap::SetInitAndTerm(unsigned int const seq_length)
{
	// |1|2|3|1|2|3|
	//	A T G          frame3, pos 0
	// 	  A T G        frame1, pos 1
	// 	    A T G      frame2, pos 2

	// initiation non-codig state

	data[0].status = iniNon;
	data[0].pos = 0;
	data[0].end_pos = 0;
	data[0].is_good = true;

	// initiation codig states: 3 frames direct and 3 frames reverse

	data[1].status = dirStart + frame3 + iniCod;
	data[2].status = revEnd   + frame3 + iniCod;
	data[3].status = dirStart + frame1 + iniCod;
	data[4].status = revEnd   + frame1 + iniCod;
	data[5].status = dirStart + frame2 + iniCod;
	data[6].status = revEnd   + frame2 + iniCod;

	data[1].pos = 0;
	data[2].pos = 0;
	data[3].pos = 1;
	data[4].pos = 1;
	data[5].pos = 2;
	data[6].pos = 2;

	//    if length%3 == 0
	//    |1|2|3|1|2|3|1|2|3|1|2|3|
	//       A T G . . . T A A
	//         A T G . . . T A A
	//           A T G . . . T A A
	
	//    if length%3 == 1
	//    |1|2|3|1|2|3|1|2|3|1|2|3|1|
	//         A T G . . . T A A
	//           A T G . . . T A A
	//             A T G . . . T A A

	//    if length%3 == 2
	//    |1|2|3|1|2|3|1|2|3|1|2|3|1|2|
	//           A T G . . . T A A
	//             A T G . . . T A A
	//               A T G . . . T A A


	unsigned int f_x = 0;
	unsigned int f_y = 0;
	unsigned int f_z = 0;

	switch ( seq_length % 3 )
	{
	case 0:
		f_x = frame1;
		f_y = frame2;
		f_z = frame3;
		break;
	case 1:
		f_x = frame2;
		f_y = frame3;
		f_z = frame1;
		break;
	case 2:
		f_x = frame3;
		f_y = frame1;
		f_z = frame2;
		break;
	}

	// termination codig states: 3 frames direct and 3 frames reverse
	
	data[data.size() - 1 - 6].status = revStart + f_x + terCod;
	data[data.size() - 1 - 5].status = dirEnd   + f_x + terCod;
	data[data.size() - 1 - 4].status = revStart + f_y + terCod;
	data[data.size() - 1 - 3].status = dirEnd   + f_y + terCod;
	data[data.size() - 1 - 2].status = revStart + f_z + terCod;
	data[data.size() - 1 - 1].status = dirEnd   + f_z + terCod;

	data[data.size() - 1 - 6].pos = seq_length - 3;
	data[data.size() - 1 - 5].pos = seq_length - 3;
	data[data.size() - 1 - 4].pos = seq_length - 2;
	data[data.size() - 1 - 3].pos = seq_length - 2;
	data[data.size() - 1 - 2].pos = seq_length - 1;
	data[data.size() - 1 - 1].pos = seq_length - 1;

	// termination non-codig state

	data[data.size() - 1 - 0].status = terNon;
	data[data.size() - 1 - 0].pos = seq_length - 1;
	data[data.size() - 1 - 0].end_pos = seq_length - 1;
	data[data.size() - 1 - 0].is_good = true;
}
// ----------------------------------------------------
void SequenceMap::SetBaseData(std::vector<unsigned int> const & flags)
{
	unsigned int total_flags = 0;

	for( vector<unsigned int>::const_iterator itr = flags.begin(); itr != flags.end(); ++itr )
	{
		if (*itr & emptyMark)
			++total_flags;
	}

	// add 7 positions for initiation and 7 positions for termination states

	total_flags += 14;

	// prepare the vector

	data.assign( total_flags, MapValue() );

	// fill in sequence-map value for nonzero flags
	// skip initiation and termination positions 
	
	unsigned int i = 7;

	for ( unsigned int j = 0; j < flags.size(); ++j )
	{
		if ( flags[j] & emptyMark )
		{
			data[i].status = flags[j];

			// point on the first NT in left codon and on last NT in right codon 
			// Atg...taA  or  Tta...caT

			if ( data[i].status & leftOrfMark )
			{
				if (data[i].status & isGap)
				{
					if (data[i].status & revEnd)
						data[i].pos = j + 3 - 2;
					else if (data[i].status & dirStart)
						data[i].pos = j + 6 - 2;
					else
						Exit("error in logic");
				}
				else
					data[i].pos = j - 2;
			}
			else if ( data[i].status & righOrftMark)
			{
				if (data[i].status & isGap)
				{
					if (data[i].status & dirEnd)
						data[i].pos = j - 3;
					else if (data[i].status & revStart)
						data[i].pos = j - 6;
					else
						Exit("error in logic");
				}
				else
					data[i].pos = j;
			}
		
			++i;
		}
	}

	// sort map by positions

	unsigned int current;
	unsigned int next;

	current = 7;
	next = 8;

	while ( next + 7 < data.size() )
	{
		if (( data[current].status & righOrftMark )&&( data[next].status & leftOrfMark ))
		{
			if ( data[current].pos >= data[next].pos )
			{
				MapValue tmp = data[current];
				data[current] = data[next];
				data[next] = tmp;
			}
		}
		
		++current;
		++next;
	}
}
// ----------------------------------------------------
void SequenceMap::Init( std::vector<unsigned int> const & flags, unsigned int const min_gene_length )
{
	data.clear();
	predictions.clear();

	SetBaseData(flags);
	SetInitAndTerm(flags.size());

	LinkDirORF(frame1);
	LinkDirORF(frame2);
	LinkDirORF(frame3);
	MarkValidDirORFs(min_gene_length);

	LinkRevORF(frame1);
	LinkRevORF(frame2);
	LinkRevORF(frame3);
	MarkValidRevORFs(min_gene_length);
}
// ----------------------------------------------------
void SequenceMap::AddCodingEvidence( std::map<int, INTERVAL_EVI > & d, std::map<int, INTERVAL_EVI > & r )
{
	vector< MapValue >::iterator itr = data.begin();
	vector< MapValue >::iterator itr_end = data.end();

	map< int, INTERVAL_EVI > evi_itr;

	for( ; itr < itr_end; ++itr )
	{
		if ( itr->status & dirEnd )
		{
			if ( d.find( itr->end_pos ) != d.end() )
				itr->IsInforced = true;
		}
		else if ( itr->status & revEnd )
		{
			if ( r.find( itr->end_pos ) != r.end() )
				itr->IsInforced = true;
		}
	}
}
// ----------------------------------------------------
void SequenceMap::LinkDirORF(unsigned int const frame)
{
	unsigned int L = -1;
	bool L_is_start = false;
	bool L_is_stop  = false;

	unsigned int R = -1;
	bool R_is_start = false;
	bool R_is_stop  = false;

	// scan all the positions from Left to Right

	for( unsigned int i = 0; i < data.size(); ++i )
	{
		if ( !(data[i].status & frame) )
			continue;

		if (data[i].status & dirStart)
		{
			R = i;
			R_is_start = true;
			R_is_stop  = false;

			if ( L_is_start )
			{
				data[L].next = &data[R];
				data[R].prev = &data[L];
			}

			L = R;
			L_is_start = R_is_start;
			L_is_stop  = R_is_stop;
		}
		else if (data[i].status & dirEnd)
		{
			R = i;
			R_is_start = false;
			R_is_stop  = true;

			if ( L_is_start )
			{
				data[L].next = &data[R];
				data[R].prev = &data[L];
			}

			L = R;
			L_is_start = R_is_start;
			L_is_stop  = R_is_stop;
		}
	}
}
// ----------------------------------------------------
void SequenceMap::LinkRevORF(unsigned int const frame)
{
	unsigned int L = -1;
	bool L_is_start = false;
	bool L_is_stop  = false;

	unsigned int R = -1;
	bool R_is_start = false;
	bool R_is_stop  = false;

	// scan all the positions from Left to Right

	for( unsigned int i = 0; i < data.size(); ++i )
	{
		if ( !(data[i].status & frame) )
			continue;

		if (data[i].status & revStart)
		{
			R = i;
			R_is_start = true;
			R_is_stop  = false;

			if ( L_is_start || L_is_stop )
			{
				data[L].next = &data[R];
				data[R].prev = &data[L];
			}

			L = R;
			L_is_start = R_is_start;
			L_is_stop  = R_is_stop;
		}
		else if (data[i].status & revEnd)
		{
			R = i;
			R_is_start = false;
			R_is_stop  = true;

			L = R;
			L_is_start = R_is_start;
			L_is_stop  = R_is_stop;
		}
	}
}
// ----------------------------------------------------
void SequenceMap::MarkValidDirORFs(unsigned int const min_length)
{
	// scan from Right to Left 

	for( vector< MapValue >::reverse_iterator itr = data.rbegin(); itr != data.rend(); ++itr )
	{
		if ( ! (itr->status & dirEnd ) )
			continue;

		itr->end_pos = itr->pos;

		MapValue* current = itr->prev;

		while( current )
		{
			current->end_pos = itr->pos;

			if ( itr->pos - current->pos + 1 >= min_length )
			{
				// found valid linked Start	
				current->is_good = true;
				itr->is_good = true;
				current = current->prev;
			}
			else
			{
				// too short: skip it in the main chain
				//    prev<-current->next
				//       C <-> B <-> A
				//       C <-------> A

				if ( current->prev )
				{
					current->next->prev = current->prev;
					current->prev->next = current->next;

					current = current->prev;
				}
				else
					current = 0;
			}
		}
	}
}
// ----------------------------------------------------
void SequenceMap::MarkValidRevORFs(unsigned int const min_length)
{
	// scan from Left to Right

	for( std::vector< MapValue >::iterator itr = data.begin(); itr != data.end(); ++itr )
	{
		if ( ! (itr->status & revEnd ) )
			continue;

		itr->end_pos = itr->pos;

		MapValue* current = itr->next;

		while( current )
		{
			current->end_pos = itr->pos;

			if ( current->pos - itr->pos + 1 >= min_length )
			{
				// found valid linked Start	
				current->is_good = true;
				itr->is_good = true;
				current = current->next;
			}
			else
			{
				// too short: skip it in the main chain
				//    prev<-current->next
				//       C <-> B <-> A
				//       C <-------> A

				if ( current->next )
				{
					current->next->prev = current->prev;
					current->prev->next = current->next;

					current = current->next;
				}
				else
					current = 0;
			}
		}
	}
}
// ----------------------------------------------------
void SequenceMap::CalcGC( Data & d )
{
	for( vector< MapValue >::iterator itr = data.begin(); itr != data.end(); ++itr )
	{
// GC is used only for is_good 
// for debug - calc for all orfs
//		if ( itr->is_good )
//		{
			if ( itr->status & dirStart )
			{
				itr->gc = d.GCasINT( itr->pos, itr->end_pos );
			}
			else if ( itr->status & revStart )
			{
				itr->gc = d.GCasINT( itr->end_pos, itr->pos );
			}
//		}
	}
}
// ----------------------------------------------------
void SequenceMap::CalcStartsGC(std::vector< Model* > & mod, std::vector<unsigned char> & nt)
{
	enum START_MODEL_TYPE { NON, SC, RBS, PROMOTER, RBS_SC, PROMOTER_SC, RBS_PROMOTER, RBS_PROMOTER_SC };
	START_MODEL_TYPE selector = NON;

	double x = 0;

	for (vector< MapValue >::iterator itr = data.begin(); itr != data.end(); ++itr)
	{
		if ((!itr->is_good) || (itr->status & iniCod) || (itr->status & terCod) || (itr->status & isGap))
			continue;

		if (itr->status & startMark)
		{
			x = mod[itr->gc]->StartContentRBS.Get(nt, itr->pos, itr->status);
			itr->logodd_start = x;
		}
	}
}
// ----------------------------------------------------
void SequenceMap::CalcStarts( Model* mod, std::vector<unsigned char> & nt )
{
	enum START_MODEL_TYPE  { NON, SC, RBS, PROMOTER, RBS_SC, PROMOTER_SC, RBS_PROMOTER, RBS_PROMOTER_SC }; 
	START_MODEL_TYPE selector = NON;

	if ( mod->RBS.is_valid && mod->StartContentRBS.is_valid && mod->Promoter.is_valid && mod->StartContentPromoter.is_valid && mod->StartContent.is_valid )
	{
		selector = RBS_PROMOTER_SC;
	}
	else if ( mod->RBS.is_valid && mod->StartContentRBS.is_valid && mod->Promoter.is_valid && mod->StartContentPromoter.is_valid )
	{
		selector = RBS_PROMOTER;
	}
	else if ( mod->Promoter.is_valid && mod->StartContentPromoter.is_valid && mod->StartContent.is_valid )
	{
		selector = PROMOTER_SC;
	}
	else if ( mod->RBS.is_valid && mod->StartContentRBS.is_valid && mod->StartContent.is_valid )
	{
		selector = RBS_SC;
	}
	else if ( mod->Promoter.is_valid && mod->StartContentPromoter.is_valid )
	{
		selector = PROMOTER;
	}
	else if ( mod->RBS.is_valid && mod->StartContentRBS.is_valid )
	{
		selector = RBS;
	}
	else if ( mod->StartContent.is_valid )
	{
		selector = SC;
	}
	else if ( !mod->RBS.is_valid && !mod->StartContentRBS.is_valid && !mod->Promoter.is_valid && !mod->StartContentPromoter.is_valid && !mod->StartContent.is_valid )
	{
		selector = NON;
	}
	else
	{
		cerr << "error, unexpected combination of start model parameters detected" << endl;
		exit(1);
	}

	double x = 0;

	for( vector< MapValue >::iterator itr = data.begin(); itr != data.end(); ++itr )
	{
		if (( !itr->is_good )||( itr->status & iniCod )||( itr->status & terCod )||(itr->status & isGap))
			continue;

		if ( itr->status & startMark )
		{
			switch (selector)
			{
				case NON:
					break;

				case SC:
					{
						itr->logodd_start = mod->StartContent.Get( nt, itr->pos, itr->status );
						itr->stype = 3;
					}
					break;

				case RBS:
					{
						x = mod->RBS.GetWithDur( nt, itr->pos, itr->status );
						itr->logodd_start = x;
						itr->logodd_RBS = x;
						x = mod->StartContentRBS.Get( nt, itr->pos, itr->status );
						itr->logodd_start += x;
						itr->logodd_RBS_SC = x;
						itr->stype = 1;
					}
					break;

				case PROMOTER:
					{
						x = mod->Promoter.GetWithDur( nt, itr->pos, itr->status );
						itr->logodd_start = x;
						itr->logodd_Promoter = x;
						x = mod->StartContentPromoter.Get( nt, itr->pos, itr->status );
						itr->logodd_start += x;
						itr->logodd_Promoter_SC = x;
						itr->stype = 2;
					}
					break;

				case RBS_SC:
					{
						double score_rbs  = mod->RBS.GetWithDur( nt, itr->pos, itr->status );
						       score_rbs += mod->StartContentRBS.Get( nt, itr->pos, itr->status );
						double score_sc  = mod->StartContent.Get( nt, itr->pos, itr->status );
						if ( score_rbs > score_sc )
						{
							itr->logodd_start = score_rbs;
							itr->stype = 1;
						}
						else
						{
							itr->logodd_start = score_sc;
							itr->stype = 3;
						}
					}
					break;

				case PROMOTER_SC:
					{
						double score_pro  = mod->Promoter.GetWithDur( nt, itr->pos, itr->status );
						       score_pro += mod->StartContentPromoter.Get( nt, itr->pos, itr->status );
						double score_sc  = mod->StartContent.Get( nt, itr->pos, itr->status );
						if ( score_pro > score_sc )
						{
							itr->logodd_start = score_pro;
							itr->stype = 2;
						}
						else
						{
							itr->logodd_start = score_sc;
							itr->stype = 3;
						}
					}
					break;

				case RBS_PROMOTER:
					{
						x = mod->RBS.GetWithDur( nt, itr->pos, itr->status );
						double score_rbs = x;
						itr->logodd_RBS = x;
						x = mod->StartContentRBS.Get( nt, itr->pos, itr->status );
						score_rbs += x;
						itr->logodd_RBS_SC = x;

						x = mod->Promoter.GetWithDur( nt, itr->pos, itr->status );
						double score_pro = x;
						itr->logodd_Promoter = x;
						x = mod->StartContentPromoter.Get( nt, itr->pos, itr->status );
						score_pro += x;
						itr->logodd_Promoter_SC = x;

						if ( score_rbs > score_pro )
						{
							itr->logodd_start = score_rbs;
							itr->stype = 1;
						}
						else
						{
							itr->logodd_start = score_pro;
							itr->stype = 2;
						}
					}
					break;

				case RBS_PROMOTER_SC:
				{
						double score_rbs  = mod->RBS.GetWithDur( nt, itr->pos, itr->status );
						       score_rbs += mod->StartContentRBS.Get( nt, itr->pos, itr->status );

						double score_pro  = mod->Promoter.GetWithDur( nt, itr->pos, itr->status );
						       score_pro += mod->StartContentPromoter.Get( nt, itr->pos, itr->status );

						double score_sc  = mod->StartContent.Get( nt, itr->pos, itr->status );

						if ( score_rbs > score_pro && score_rbs > score_sc )
						{
							itr->logodd_start = score_rbs;
							itr->stype = 1;
						}
						else if ( score_pro > score_rbs && score_pro > score_sc )
						{
							itr->logodd_start = score_pro;
							itr->stype = 2;
						}
						else if ( score_sc > score_rbs && score_sc > score_pro )
						{
							itr->logodd_start = score_sc;
							itr->stype = 3;
						}
						else
						{
							cerr << "error, check code" << endl;
							exit(1);
						}
					}
					break;
			}
		}
	}
}
// ----------------------------------------------------
void SequenceMap::CalcLogP( std::vector< Model* > & mod, std::vector<unsigned char> & nt, char const gtype )
{
	vector< MapValue >::iterator  itr     = data.begin();
	vector< MapValue >::iterator  itr_end = data.end();

	for( ; itr != itr_end; ++itr )
	{
		if ( !itr->is_good )               continue;
		if ( !(itr->status & startMark) )  continue;

		vector<double> * logP_1     = & (mod[ itr->gc ]->logP_1_n);
		vector<double> * logP_2     = & (mod[ itr->gc ]->logP_2_n);
		vector<double> * logP_3     = & (mod[ itr->gc ]->logP_3_n);
		vector<double> * logP_1_abs = & (mod[ itr->gc ]->logP_1_abs_n);

		unsigned int order =  mod[ itr->gc ]->order_cod;

		Index idx(order);
		unsigned int current_index;

		// do calc here

		if ( itr->status & dirStart )
		{
			unsigned int current_pos = itr->pos;
				
			if ((itr->status & iniCod)||(itr->status & isGap))
				;
			else
				current_pos = current_pos + mod[ itr->gc ]->ORF_start_marging;

			for( unsigned int i = 0; i <= order; ++i )
			{
				current_index = idx.GetWithN( nt[ current_pos ] );
				++current_pos;
			}

			double logP = logP_1_abs->operator[](current_index);

			current_index = idx.GetWithN( nt[current_pos] );
			logP += logP_2->operator[](current_index);
			++current_pos;

			current_index = idx.GetWithN( nt[current_pos] );
			logP += logP_3->operator[](current_index);
			++current_pos;

			unsigned int end_pos = itr->end_pos - 2;

			while(1)
			{
				if ( current_pos < end_pos )
				{
					current_index = idx.GetWithN( nt[current_pos] );
					logP += logP_1->operator[](current_index);
					++current_pos;
				}
				else break;

				if ( current_pos < end_pos )
				{
					current_index = idx.GetWithN( nt[current_pos] );
					logP += logP_2->operator[](current_index);
					++current_pos;
				}
				else break;

				if ( current_pos < end_pos )
				{
					current_index = idx.GetWithN( nt[current_pos] );
					logP += logP_3->operator[](current_index);
					++current_pos;
				}
				else break;
			}

			logP += mod[ itr->gc ]->toModel;

			if (  itr->logP < logP  )
			{
				itr->logP = logP;
				itr->logP_type = gtype;
			}
		}
		else if ( itr->status & revStart )
		{
			unsigned int current_pos = itr->pos;
				
			if ((itr->status & iniCod) || (itr->status & isGap))
				;
			else
				current_pos = current_pos - mod[ itr->gc ]->ORF_start_marging;

			for ( unsigned int i = 0; i <= order; ++i )
			{
				current_index = idx.GetRCwithN( nt[ current_pos] );
				--current_pos;
			}

			double logP = logP_1_abs->operator[](current_index);

			current_index = idx.GetRCwithN( nt[current_pos] );
			logP += logP_2->operator[](current_index);
			--current_pos;

			current_index = idx.GetRCwithN( nt[current_pos] );
			logP += logP_3->operator[](current_index);
			--current_pos;

			unsigned int end_pos = itr->end_pos + 2;

			while (1)
			{
				if ( current_pos > end_pos )
				{
					current_index = idx.GetRCwithN( nt[current_pos] );
					logP += logP_1->operator[](current_index);
					--current_pos;
				}
				else break;

				if ( current_pos > end_pos )
				{
					current_index = idx.GetRCwithN( nt[current_pos] );
					logP += logP_2->operator[](current_index);
					--current_pos;
				}
				else break;

				if ( current_pos > end_pos )
				{
					current_index = idx.GetRCwithN( nt[current_pos] );
					logP += logP_3->operator[](current_index);
					--current_pos;
				}
				else break;
			}

			logP += mod[ itr->gc ]->toModel;

			if (  itr->logP < logP  )
			{
				itr->logP = logP;
				itr->logP_type = gtype;
			}
		}
	}
}
// ----------------------------------------------------
void SequenceMap::CalcLogodd( std::vector< Model* > & mod, std::vector<unsigned char> & nt, char const gtype )
{
	for( vector< MapValue >::iterator itr = data.begin(); itr != data.end(); ++itr )
	{
		if ( !itr->is_good )               continue;
		if ( !(itr->status & startMark) )  continue;

		vector<double> * logodd_1     = & (mod[ itr->gc ]->logodd_1);
		vector<double> * logodd_2     = & (mod[ itr->gc ]->logodd_2);
		vector<double> * logodd_3     = & (mod[ itr->gc ]->logodd_3);
		vector<double> * logodd_1_abs = & (mod[ itr->gc ]->logodd_1_abs);

		unsigned int order =  mod[ itr->gc ]->order_cod; 

		Index idx(order);
		unsigned int current_index;

		// temporary check point 

		if ( mod[ itr->gc ]->ORF_start_marging % 3 )
		{
			Exit("error, only 'ORF_start_marging % 3 == 0 ' is supported ");
		}

		if ( itr->status & dirStart )
		{
			unsigned int current_pos = itr->pos;
				
			if ((itr->status & iniCod) || (itr->status & isGap))
				;
			else
				current_pos = current_pos + mod[ itr->gc ]->ORF_start_marging;

			for( unsigned int i = 0; i <= order; ++i )
			{
				current_index = idx.GetWithN( nt[ current_pos ] );
				++current_pos;
			}

			double logodd = logodd_1_abs->operator[](current_index);

			current_index = idx.GetWithN( nt[current_pos] );
			logodd += logodd_2->operator[](current_index);
			++current_pos;

			current_index = idx.GetWithN( nt[current_pos] );
			logodd += logodd_3->operator[](current_index);
			++current_pos;

			unsigned int end_pos = itr->end_pos - 2;

			while(1)
			{
				if ( current_pos < end_pos )
				{
					current_index = idx.GetWithN( nt[current_pos] );
					logodd += logodd_1->operator[](current_index);
					++current_pos;
				}
				else break;

				if ( current_pos < end_pos )
				{
					current_index = idx.GetWithN( nt[current_pos] );
					logodd += logodd_2->operator[](current_index);
					++current_pos;
				}
				else break;

				if ( current_pos < end_pos )
				{
					current_index = idx.GetWithN( nt[current_pos] );
					logodd += logodd_3->operator[](current_index);
					++current_pos;
				}
				else break;
			}

			MapValue* ptr_to_end = itr->next;

			while( ptr_to_end->next )
				ptr_to_end = ptr_to_end->next;

			double start_stop = 0;

			if      ( ptr_to_end->status & isTAA )  start_stop += mod[ itr->gc ]->logodd_TAA;
			else if ( ptr_to_end->status & isTAG )  start_stop += mod[ itr->gc ]->logodd_TAG;
			else if ( ptr_to_end->status & isTGA )  start_stop += mod[ itr->gc ]->logodd_TGA;

			if      ( itr->status & isATG )  start_stop += mod[ itr->gc ]->logodd_ATG;
			else if ( itr->status & isGTG )  start_stop += mod[ itr->gc ]->logodd_GTG;
			else if ( itr->status & isTTG )  start_stop += mod[ itr->gc ]->logodd_TTG;

			double dur = mod[ itr->gc ]->GetDurationForORF( itr->end_pos - itr->pos + 1 );

			if ( itr->status & incomplete_ORF || ptr_to_end->status & incomplete_ORF )
			{
				dur = mod[ itr->gc ]->GetIncompleteDurationForORF( itr->end_pos - itr->pos + 1 );
			}

			double current_total = logodd + start_stop + dur + itr->logodd_start + mod[itr->gc]->toModel;

			if ( ptr_to_end->IsInforced )
			{
				itr->IsInforced = true;
				current_total += 100;  // change to score
			}

			if ( current_total >= 0 )
				itr->gtype += ( gtype << 4 );

			if ( itr->logP_type == gtype )
//			if (  itr->prob < current_total  )
			{
				itr->logodd_orf = logodd;
				itr->logodd_dur = dur;
				itr->prob = current_total;
				itr->gtype &= 240;
				itr->gtype += gtype;
			}

			// tmp code, development step
			
			{
				if ( gtype == NATIVE_TYPE )
				{
					itr->logodd_total_native = current_total;
					itr->logodd_orf_native   = logodd;
					itr->logodd_dur_native   = dur;
					itr->logodd_tr_native    = mod[itr->gc]->toModel;
				}
				else
				{
					if ( current_total > itr->logodd_total_atypical )
					{
						itr->logodd_total_atypical = current_total;
						itr->logodd_orf_atypical   = logodd;
						itr->logodd_dur_atypical   = dur;
						itr->logodd_tr_atypical    = mod[itr->gc]->toModel;

					}
				}

				itr->logodd_start_stop = start_stop;
				itr->logodd_start = itr->logodd_start;
			}

			// end of tmp code, development step
		}
		else if ( itr->status & revStart )
		{
			unsigned int current_pos = itr->pos;
				
			if ((itr->status & iniCod) || (itr->status & isGap))
				;
			else
				current_pos = current_pos - mod[ itr->gc ]->ORF_start_marging;

			for ( unsigned int i = 0; i <= order; ++i )
			{
				current_index = idx.GetRCwithN( nt[ current_pos] );
				--current_pos;
			}

			double logodd = logodd_1_abs->operator[](current_index);

			current_index = idx.GetRCwithN( nt[current_pos] );
			logodd += logodd_2->operator[](current_index);
			--current_pos;

			current_index = idx.GetRCwithN( nt[current_pos] );
			logodd += logodd_3->operator[](current_index);
			--current_pos;

			unsigned int end_pos = itr->end_pos + 2;

			while (1)
			{
				if ( current_pos > end_pos )
				{
					current_index = idx.GetRCwithN( nt[current_pos] );
					logodd += logodd_1->operator[](current_index);
					--current_pos;
				}
				else break;

				if ( current_pos > end_pos )
				{
					current_index = idx.GetRCwithN( nt[current_pos] );
					logodd += logodd_2->operator[](current_index);
					--current_pos;
				}
				else break;

				if ( current_pos > end_pos )
				{
					current_index = idx.GetRCwithN( nt[current_pos] );
					logodd += logodd_3->operator[](current_index);
					--current_pos;
				}
				else break;
			}

			MapValue* ptr_to_end = itr->prev;

			while( ptr_to_end->prev )
				ptr_to_end = ptr_to_end->prev;

			double start_stop = 0;

			if      ( ptr_to_end->status & isTAA )  start_stop += mod[ itr->gc ]->logodd_TAA;
			else if ( ptr_to_end->status & isTAG )  start_stop += mod[ itr->gc ]->logodd_TAG;
			else if ( ptr_to_end->status & isTGA )  start_stop += mod[ itr->gc ]->logodd_TGA;

			if      ( itr->status & isATG )  start_stop += mod[ itr->gc ]->logodd_ATG;
			else if ( itr->status & isGTG )  start_stop += mod[ itr->gc ]->logodd_GTG;
			else if ( itr->status & isTTG )  start_stop += mod[ itr->gc ]->logodd_TTG;

			double dur = mod[ itr->gc ]->GetDurationForORF( itr->pos - itr->end_pos + 1 );

			if ( itr->status & incomplete_ORF || ptr_to_end->status & incomplete_ORF )
			{
				dur = mod[ itr->gc ]->GetIncompleteDurationForORF( itr->pos - itr->end_pos + 1);
			}

			double current_total = logodd + start_stop + dur + itr->logodd_start + mod[itr->gc]->toModel;

			if ( ptr_to_end->IsInforced )
			{
				itr->IsInforced = true;
				current_total += 100;  // change to score
			}

			if ( current_total >= 0 )
				itr->gtype += ( gtype << 4 );

			if ( itr->logP_type == gtype )
//			if (  itr->prob < current_total )
			{
				itr->logodd_orf = logodd;
				itr->logodd_dur = dur;
				itr->prob = current_total;
				itr->gtype &= 240;
				itr->gtype += gtype;
			}

			// tmp code, development step
			
			{
				if ( gtype == NATIVE_TYPE )
				{
					itr->logodd_total_native = current_total;
					itr->logodd_orf_native   = logodd;
					itr->logodd_dur_native   = dur;
					itr->logodd_tr_native    = mod[itr->gc]->toModel;
				}
				else
				{
					if ( current_total > itr->logodd_total_atypical )
					{
						itr->logodd_total_atypical = current_total;
						itr->logodd_orf_atypical   = logodd;
						itr->logodd_dur_atypical   = dur;
						itr->logodd_tr_atypical    = mod[itr->gc]->toModel;
					}
				}

				itr->logodd_start_stop = start_stop;
				itr->logodd_start = itr->logodd_start;
			}

			// end of tmp code, development step
		}
	}
}
// ----------------------------------------------------
void SequenceMap::ExcludeNegativeScoredDirORFs(double delta)
{
	for( std::vector< MapValue >::iterator itr = data.begin(); itr != data.end(); ++itr )
	{
		if ( itr->is_good && ( itr->status & dirEnd ) )
		{
			for( MapValue* current = itr->prev; current != 0; current = current->prev )
			{
				if ( current->status & dirStart )
				{
					current->prob += delta;

					if ( current->prob < 0 )
					{
						current->is_good = false;
					
						if ( current->prev )
						{
							current->next->prev = current->prev;
							current->prev->next = current->next;
						}
						else
						{
							current->next->prev = 0;
						}
					}
				}
				else
				{
					cerr << "error in 'SequenceMap::ExcludenegativeScoredDirORFs()'" << endl;
					exit(1);
				}
			}

			if ( ! itr->prev )
				itr->is_good = false;
		}
	}
}
// ----------------------------------------------------
void SequenceMap::ExcludeNegativeScoredRevORFs(double delta)
{
	for( std::vector< MapValue >::iterator itr = data.begin(); itr != data.end(); ++itr )
	{
		if ( itr->is_good && ( itr->status & revEnd ) )
		{
			for( MapValue* current = itr->next; current != 0; current = current->next )
			{
				if ( current->status & revStart )
				{
					current->prob += delta;

					if ( current->prob < 0 )
					{
						current->is_good = false;
					
						if ( current->next )
						{
							current->next->prev = current->prev;
							current->prev->next = current->next;
						}
						else
						{
							current->prev->next = 0;
						}
					}
				}
				else
				{
					cerr << "error in 'SequenceMap::ExcludenegativeScoredRevORFs()'" << endl;
					exit(1);
				}
			}

			if ( ! itr->next )
				itr->is_good = false;
		}
	}
}
// ----------------------------------------------------
void SequenceMap::FindMaxScoredORFs(void)
{
	FindMaxScoredDirORFs();
	FindMaxScoredRevORFs();
//	TestORFsDir();
//	TestORFsRev();
}
// ----------------------------------------------------
void SequenceMap::CalcBayesDirORFs(void)
{
	for (std::vector< MapValue >::iterator itr = data.begin(); itr != data.end(); ++itr)
	{
		double max = 0;

		// find max score for stop codon
		if (itr->is_good && (itr->status & revEnd))
		{
			// Find max score
			for (MapValue* current = &(*itr); current != 0; current = current->next)
			{
				if (current->status & revStart)
				{
					if (current->prob > max)
					{
						max = current->prob;
					}
				}
			}

			if (max <= 0)
			{
				continue;
			}

			double sum = 0;

			// calc norm
			for (MapValue* current = &(*itr); current != 0; current = current->next)
			{
				if (current->status & revStart)
				{ 
					if (current->prob - max > -20)
					{
						sum += exp(current->prob - max);
					}
				}
			}

			// calc Bayes
			for (MapValue* current = &(*itr); current != 0; current = current->next)
			{
				if (current->status & revStart)
				{
					if (current->prob - max > -20)
					{
						current->bayes = exp(current->prob - max)/sum;
					}
					else
					{
						current->bayes =  0;
					}
				}
			}
		}
	}
}
// ----------------------------------------------------
void SequenceMap::CalcBayesRevORFs(void)
{
	for (std::vector< MapValue >::iterator itr = data.begin(); itr != data.end(); ++itr)
	{
		double max = 0;

		// find max score for stop codon
		if (itr->is_good && (itr->status & dirEnd))
		{
			// Find max score
			for (MapValue* current = &(*itr); current != 0; current = current->prev)
			{
				if (current->status & dirStart)
				{
					if (current->prob > max)
					{
						max = current->prob;
					}
				}
			}

			if (max <= 0)
			{
				continue;
			}

			double sum = 0;

			// calc norm
			for (MapValue* current = &(*itr); current != 0; current = current->prev)
			{
				if (current->status & dirStart)
				{
					if (current->prob - max > -20)
					{
						sum += exp(current->prob - max);
					}
				}
			}

			// calc Bayes
			for (MapValue* current = &(*itr); current != 0; current = current->prev)
			{
				if (current->status & dirStart)
				{
					if (current->prob - max > -20)
					{
						current->bayes = exp(current->prob - max) / sum;
					}
					else
					{
						current->bayes = 0;
					}
				}
			}
		}
	}
}
// ----------------------------------------------------
void SequenceMap::FindMaxScoredDirORFs(void)
{
	for( std::vector< MapValue >::iterator itr = data.begin(); itr != data.end(); ++itr )
	{
		if ( itr->is_good &&( itr->status & dirStart ))
		{
			double max = 0;

			for( MapValue* current = &(*itr); current != 0; current = current->next )
			{
				if ( current->status & dirStart )
				{
					if ( current->prob > max )
					{
						max = current->prob;

						if ( current->prev )
						{
							current->prev->is_good = false;
							current->prev = 0;
						}
					}
					else
					{
						current->is_good = false;

						if ( current->prev )
						{
							current->prev->next = current->next;
							current->next->prev = current->prev;
						}
					}
				}
				else if ( current->status & dirEnd )
				{
					if ( max == 0 )
						current->is_good = false;
				}
			}
		}
	}
}
// ----------------------------------------------------
void SequenceMap::FindMaxScoredRevORFs(void)
{
	for( std::vector< MapValue >::reverse_iterator itr = data.rbegin(); itr != data.rend(); ++itr )
	{
		if ( itr->is_good &&( itr->status & revStart ))
		{
			double max = 0;

			for( MapValue* current = &(*itr); current != 0; current = current->prev )
			{
				if ( current->status & revStart )
				{
					if ( current->prob > max )
					{
						max = current->prob;

						if ( current->next )
						{
							current->next->is_good = false;
							current->next = 0;
						}
					}
					else
					{
						current->is_good = false;

						if ( current->next )
						{
							current->next->prev = current->prev;
							current->prev->next = current->next;
						}
					}
				}
				else if ( current->status & revEnd )
				{
					if ( max == 0 )
						current->is_good = false;
				}
			}
		}
	}
}
// ----------------------------------------------------
std::string SequenceMap::PrintMap(void)
{
	stringstream s;
	int count_status = 0; 

	s << "## map information current record" << endl;
	s << "#1 index of point in vector of sequence map" << endl;
	s << "#2 TRUE/FALSE status of point in sequence map" << endl;
	s << "#3 ORF G+C for start codon location" << endl;
	s << "#4 position on sequence" << endl;
	s << "#5 position STOP codon for current START codon" << endl;
	s << "#6 length of ORF for STARTS" << endl;
	s << "#7 TYPE of flag" << endl;
	s << "#8 position on sequence of prev" << endl;
	s << "#9 position on sequence of next" << endl;
	s << "#10 total loggodd best" << endl;
	s << "#11 total loggodd native" << endl;
	s << "#12 total loggodd atypical" << endl;
	s << "#11 loggodd ORF native" << endl;
	s << "#12 loggodd ORF atypical" << endl;
	s << "#13 loggodd duration native" << endl;
	s << "#14 loggodd duration atypical" << endl;
	s << "#15 loggodd TR native" << endl;
	s << "#16 loggodd TR atypical" << endl;
	s << "#17 loggodd START total" << endl;
	s << "#18 loggodd start + stop" << endl;
	s << "#19 loggodd RBS" << endl;
	s << "#20 loggodd RBS SC" << endl;
	s << "#21 loggodd Promoter" << endl;
	s << "#22 loggodd Promoter SC" << endl;

	for ( unsigned int i = 0; i < data.size(); ++i )
	{	
		if ( !data[i].is_good ) continue;

		++count_status;

		s << i << " " << data[i].is_good << " ";

		if ( data[i].status & startMark )
			s << data[i].gc << " ";
		else
			s << 0 << " ";

		if ( data[i].end_pos > data[i].pos )
		{
			s << data[i].pos << " " << data[i].end_pos << " " << data[i].end_pos - data[i].pos + 1 << " ";
		}
		else if ( data[i].end_pos < data[i].pos )
		{
			s << data[i].pos << " " << data[i].end_pos << " " << data[i].pos - data[i].end_pos + 1 << " ";
		}
		else
			s << data[i].pos << " 0 0 ";

		if ( data[i].status & dirEnd )   s << "D-end_";
		if ( data[i].status & revEnd )   s << "R-end_";
		if ( data[i].status & dirStart ) s << "D-sta_";
		if ( data[i].status & revStart ) s << "R-sta_";
		if ( data[i].status & frame1 )   s << "fr1_";
		if ( data[i].status & frame2 )   s << "fr2_";
		if ( data[i].status & frame3 )   s << "fr3_";
		if ( data[i].status & isATG )    s << "atg_";
		if ( data[i].status & isGTG )    s << "gtg_";
		if ( data[i].status & isTTG )    s << "ttg_";
		if ( data[i].status & isTAA )    s << "taa_";
		if ( data[i].status & isTAG )    s << "tag_";
		if ( data[i].status & isTGA )    s << "tga_";
		if ( data[i].status & iniCod )   s << "ini-C_";
		if ( data[i].status & iniNon )   s << "ini-N_";
		if ( data[i].status & terCod )   s << "ter-C_";
		if ( data[i].status & terNon )   s << "ter-N_";

		if ( data[i].prev )
			s << " " << data[i].prev->pos << " ";
		else
			s << " " << 0 << " ";
		
		if ( data[i].next )
			s << data[i].next->pos << " ";
		else
			s << "0" << " ";

		s << "; ";

		s << data[i].prob << " " << data[i].logodd_total_native << " " << data[i].logodd_total_atypical << " ";
		s << data[i].logodd_orf_native << " " << data[i].logodd_orf_atypical << " ";
		s << data[i].logodd_dur_native << " " << data[i].logodd_dur_atypical << " ";
		s << data[i].logodd_tr_native << " " << data[i].logodd_tr_atypical << " ";

		s << data[i].logodd_start << " " << data[i].logodd_start_stop << " ";
		s << data[i].logodd_RBS << " " << data[i].logodd_RBS_SC << " ";
		s << data[i].logodd_Promoter << " " << data[i].logodd_Promoter_SC << " ";

		s << "; ";

		s << data[i].vprob << " ";

		if ( data[i].path )
			s << data[i].path->pos;
		else
			s << "0" ;

		if ( data[i].gtype & NATIVE_TYPE )
			s << " native";
		if ( data[i].gtype & ATYPICAL_TYPE_1 )
			s << " bac";
		if ( data[i].gtype & ATYPICAL_TYPE_2 )
			s << " arc";

		s << endl;
	}

	s << "total non zero map positions: " << count_status << endl;

	return s.str();
}
// ----------------------------------------------------
void SequenceMap::DP_initialization_of_noncoding_state(vector< MapValue >::iterator itr)
{
	if ( !(itr->status & iniNon) )
	{
		cout << "error 'SequenceMap::DP_initialization_of_noncoding_state()'" << endl;
		exit(1);
	}

	itr->path = 0;
	itr->vprob = 0;
}
// ----------------------------------------------------
void SequenceMap::DP_initialization_of_coding_state(vector< MapValue >::iterator itr)
{
	if ( !(itr->status & iniCod) )
	{
		cout << "error 'SequenceMap::DP_initialization_of_coding_state()'" << endl;
		exit(1);
	}

	itr->path = 0;
	itr->vprob = 0;
}
// ----------------------------------------------------
void SequenceMap::DP_ORF(vector< MapValue >::iterator itr)
{
	if ( itr->status & dirEnd )
	{
		MapValue* other_end = itr->prev;

		while( other_end )
		{
			SelectBestCodDir( other_end, itr );

			vector< MapValue >::iterator ov = itr;

			do
			{
				--ov;

				if ( ov->is_good &&( ov->status & righOrftMark )&& ov->path && ov->path < other_end)
				{
					if ((itr->status & terCod) && (ov->status & terCod))
						continue;

					if ((itr->status & isGap) && (ov->status & isGap))
						continue;

					if ((other_end->status & iniCod) && (ov->path->status & iniCod))
						continue;

					if ((other_end->status & isGap) && (ov->path->status & isGap))
						continue;

					SelectBestOverlapCodDir( other_end, ov, itr );
				}
			} 
			while( &(*ov) > other_end );

			other_end = other_end->prev;
		}
	}
	else if ( itr->status & revStart )
	{
		MapValue* other_end = itr->prev;

		while( other_end->prev )
		{
			other_end = other_end->prev;
		}

		SelectBestCodRev( other_end, itr );

		vector< MapValue >::iterator ov = itr;

		do
		{
			--ov;
			
			if ( ov->is_good &&( ov->status & righOrftMark )&& ov->path && ov->path < other_end)
			{
				if ((itr->status & terCod) && (ov->status & terCod))
					continue;

				if ((itr->status & isGap) && (ov->status & isGap))
					continue;

				if ((other_end->status & iniCod) && (ov->path->status & iniCod))
					continue;

				if ((other_end->status & isGap) && (ov->path->status & isGap))
					continue;

				SelectBestOverlapCodRev( other_end, ov, itr );
			}
		}
		while( &(*ov) > other_end );
	}
}
// ----------------------------------------------------
void SequenceMap::DP_intergenic(std::vector< MapValue >::iterator itr)
{
	vector< MapValue >::iterator prev = itr;

	int counter = 100;

	do
	{
		--prev;
		
		if ( prev->is_good && (( prev->status & codMark )||( prev->status & iniNon )))
		{
			if (prev->status & terCod)
				continue;

			SelectBestN( prev, itr );
			--counter;
		}
	}
	while(( prev != data.begin() ) && counter );
}
// ----------------------------------------------------
void SequenceMap::Dynamic(void)
{
	vector< MapValue >::iterator current = data.begin(); 

	DP_initialization_of_noncoding_state( current );
	++current;

	while(( current != data.end() ) && ( current->status & iniCod ))
	{
		if ( current->is_good )
			DP_initialization_of_coding_state( current );

		++current;
	}

	while( current < data.end() )
	{
		if ( !current->is_good )
		{
			;
		}
		else if ( current->status & codMark )
		{
			DP_ORF(current);
		}
		else if (( current->status & noncodMark )||( current->status & terNon ))
		{
			DP_intergenic(current);
		}
		else
		{
			cout << "error 'SequenceMap::Dynamic()'" << endl;
			exit(1);
		}

		++current;
	}
}
// ----------------------------------------------------
void SequenceMap::BestPath(void)
{
	vector< MapValue >::reverse_iterator best;
	double best_end_value = LOG_ZERO;

	vector< MapValue >::reverse_iterator current = data.rbegin();

	if ( current->vprob > best_end_value )
	{
		best_end_value = current->vprob;
		best = current;
	}

	++current;

	while( current != data.rend() && (current->status & terCod) )
	{
		if ( current->vprob > best_end_value )
		{
			best_end_value = current->vprob;
			best = current;
		}

		++current;
	}

	final_logodd = best_end_value;

	// estimate size of best path
	current = best;
	MapValue* p = &(*current);

	unsigned int best_path_size = 0;
	
	while(p)
	{
		++best_path_size;
		p = p->path;
	}

	// get best path

	predictions.reserve( best_path_size );

	current = best;
	p = &(*current);

	while(p)
	{
		BestValue best_value;

		if ( p->gtype )
		{
			if (p->end_pos > p->pos)
			{
				best_value.L = p->pos + 1;
				best_value.R = p->end_pos + 1;
				best_value.strand = STRAND_TYPES::DIRECT;
				best_value.gtype = p->gtype;

				best_value.complete_L = true;
				best_value.complete_R = true;

				if (p->status & iniCod || p->status & isGap)
					best_value.complete_L = false;

				MapValue* p_end = p->next;
				while (p_end->next)
					p_end = p_end->next;
				if (p_end->status & terCod || p_end->status & isGap)
					best_value.complete_R = false;

				best_value.stype = p->stype;

				best_value.origin = p;

				predictions.push_back(best_value);
			}
			else if ( p->end_pos < p->pos )
			{
				best_value.L = p->end_pos + 1;
				best_value.R = p->pos + 1;
				best_value.strand = STRAND_TYPES::REVERSE;
				best_value.gtype = p->gtype;

				best_value.complete_L = true;
				best_value.complete_R = true;

				if (p->status & terCod || p->status & isGap)
					best_value.complete_R = false;

				if (p->path->status & iniCod || p->path->status & isGap)
					best_value.complete_L = false;

				best_value.stype = p->stype;

				best_value.origin = p;
				predictions.push_back(best_value);
			}
		}

		p = p->path;
	}
}
// ----------------------------------------------------
void SequenceMap::SelectBestN( std::vector< MapValue >::iterator p, std::vector< MapValue >::iterator c )
{
	if ( c->vprob < p->vprob )
	{
		c->vprob = p->vprob;
		c->path = &(*p);
	}
}
// ----------------------------------------------------
void SequenceMap::SelectBestCodDir( MapValue* p, std::vector< MapValue >::iterator c )
{
	// score is stored at start position

	if ( c->vprob < p->vprob + p->prob )
	{
		c->vprob = p->vprob + p->prob;
		c->path = p;
	}
}
// ----------------------------------------------------
void SequenceMap::SelectBestCodRev( MapValue* p, std::vector< MapValue >::iterator c )
{
	if ( c->vprob < p->vprob + c->prob )
	{
		c->vprob = p->vprob + c->prob;
		c->path = p;
	}
}
// ----------------------------------------------------
void SequenceMap::SelectBestOverlapCodDir(MapValue* p, std::vector< MapValue >::iterator o, std::vector< MapValue >::iterator c)
{
	double y = c->pos - p->pos + 1.0;
	double z = o->pos - o->path->pos + 1.0;
	double x = o->pos - p->pos + 1.0;

	double penalty = 0;

	if (o->status & dirORF)
	{
		// direct to direct overlap

		if (z < y)
			penalty = x*log(1 + 0.4 * overlap_penalty*(x / y + x / z));
		else
			penalty = x*log(1 + overlap_penalty*(x / y + x / z));
	}
	else if (o->status & revORF)
	{
		// direct to reverse overlap : start to start - larger penalty

		penalty = x*log(1 + 2 * overlap_penalty*(x / y + x / z));
	}
	else
	{
		cerr << "error, 'SequenceMap::SelectBestOverlapCodDir()'" << endl;
		exit(1);
	}

	if ( c->vprob < p->prob + o->vprob - penalty )
	{
		c->vprob = p->prob + o->vprob - penalty;

		c->path = p;
		p->path = &(*o);
	}
}
// ----------------------------------------------------
void SequenceMap::SelectBestOverlapCodRev( MapValue* p, std::vector< MapValue >::iterator o, std::vector< MapValue >::iterator c )
{
	double y = c->pos - p->pos + 1.0; 
	double z = o->pos - o->path->pos + 1.0;
	double x = o->pos - p->pos + 1.0;

	double penalty = 0;
	
	if ( o->status & revORF )
	{
		// reverse to reverse overlap

		if ( y < z )
			penalty = x*log( 1 + 0.4 * overlap_penalty*( x/y + x/z) );
		else
			penalty = x*log( 1 + overlap_penalty*( x/y + x/z) );
	}
	else if ( o->status & dirORF )
	{
		// reverse to direct overlap : end to end

		penalty = x*log( 1 + overlap_penalty*( x/y + x/z) );
	}
	else 
	{
		cerr << "error, 'SequenceMap::SelectBestOverlapCodRev()'" << endl;
		exit(1);
	}

	if ( c->vprob < c->prob + o->vprob - penalty )
	{
		c->vprob = c->prob + o->vprob - penalty;
		c->path = p;
		p->path = &(*o);
	}
}
// ----------------------------------------------------
void SequenceMap::ParseIsland( std::vector< MapValue >::iterator L, std::vector< MapValue >::iterator R )
{
	;
}
// ----------------------------------------------------
bool SequenceMap::FindIsland( std::vector< MapValue >::iterator & L, std::vector< MapValue >::iterator & R )
{
	vector< MapValue >::iterator current = L;

	while( current < R )
	{
		if ( current->is_good &&( current->status & noncodMark ) )
		{
			break;
		}
		else
			++current;
	}

	if ( current == R )
		return false;

	L = current;

	// pointer to iterator
	vector< MapValue >::iterator linked = data.begin() + ( current->next - &(*data.begin()) );
	
	++current;

	while( current < R )
	{
		if ( current == linked )
			break;

		if ( current->is_good )
		{
			if ( current->next && ( current->next > &(*linked) ))
			{
				linked = data.begin() + ( current->next - &(*data.begin()) );
			}
		}

		++current;
	}

	if ( current == R )
		return false;

	R = current;
	
	return true;
}
// ----------------------------------------------------
void SequenceMap::TestORFsDir(void)
{
	for( vector< MapValue >::iterator itr = data.begin(); itr != data.end(); ++itr )
	{
		if ( itr->is_good )
		{
			if ( itr->status & dirStart )
			{
				if (  itr->next->is_good && ( itr->next->status & dirEnd ) )
				{
					itr->next->prob = itr->prob;
					itr->next->end_pos = itr->pos;
					itr->next->gc = itr->gc;
				}
				else
				{
					cerr << "error with dir stop" << endl;
				}
			}

			if ( itr->status & dirEnd )
			{
				if (  itr->prev->is_good && ( itr->prev->status & dirStart ) )
					;
				else
				{
					cerr << "error with dir start" << endl;
				}
			}
		}
	}
}
// ----------------------------------------------------
void SequenceMap::TestORFsRev(void)
{
	for( vector< MapValue >::iterator itr = data.begin(); itr != data.end(); ++itr )
	{
		if ( itr->is_good )
		{
			if ( itr->status & revStart )
			{
				if (  itr->prev->is_good && ( itr->prev->status & revEnd ) )
				{
					itr->prev->prob = itr->prob;
					itr->prev->end_pos = itr->pos;
					itr->prev->gc = itr->gc;
				}
				else
				{
					cerr << "error with rev stop" << endl;
				}
			}

			if ( itr->status & revEnd )
			{
				if (  itr->next->is_good && ( itr->next->status & revStart ) )
					;
				else
				{
					cerr << "error with rev start" << endl;
				}
			}
		}
	}
}
// ----------------------------------------------------
void SequenceMap::Run(bool best_start_before_dp, double delta)
{
	AddIniTermProb( noncoding_state_initiation_and_termination_probability );

	ExcludeNegativeScoredDirORFs(delta);
	ExcludeNegativeScoredRevORFs(delta);

	if ( best_start_before_dp )
		FindMaxScoredORFs();

	CalcBayesDirORFs();
	CalcBayesRevORFs();

	Dynamic();
	BestPath();

	if (logger->debug) logger->Print( 1, PrintMap() );
}
// ----------------------------------------------------
void SequenceMap::AddSiteInfo( Model* mod, std::vector<unsigned char> & nt )
{
	for( vector< BestValue >::iterator itr = predictions.begin(); itr != predictions.end(); ++itr )
	{
		if ( itr->strand == STRAND_TYPES::DIRECT )
		{
			if (itr->stype == 1)
			{
				itr->si = mod->RBS.GetDirWithDurFullInfo( nt, itr->L - 1 );
				itr->si.stype = 1;
			}
			else if (itr->stype == 2)
			{
				itr->si = mod->Promoter.GetDirWithDurFullInfo( nt, itr->L - 1 );
				itr->si.stype = 2;
			}
			else if (itr->stype == 3)
			{
				itr->si.stype = 3;
			}
		}
		else if ( itr->strand == STRAND_TYPES::REVERSE )
		{
			if (itr->stype == 1)
			{
				itr->si = mod->RBS.GetRevWithDurFullInfo( nt, itr->R - 1 );
				itr->si.stype = 1;
			}
			else if (itr->stype == 2)
			{
				itr->si = mod->Promoter.GetRevWithDurFullInfo( nt, itr->R - 1 );
				itr->si.stype = 2;
			}
			else if (itr->stype == 3)
			{
				itr->si.stype = 3;
			}
		}
	}
}
// ----------------------------------------------------

