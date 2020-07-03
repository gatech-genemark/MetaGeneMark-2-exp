// ====================================================
// Author: Alex Lomsadze
// Copyright: Georgia Institute of Technology
// Copyright: GeneProbe Inc.
// Project: GeneMark.hmm-2 (no introns)
// File: model_2.h
// ====================================================

#ifndef MODEL_2_H
#define MODEL_2_H

#include <string>
#include <vector>

#include "logger.h"
#include "common_2.h"
#include "site_2.h"

// ----------------------------------------------------
namespace DEFVALUES
{
	double const ERROR =  -1.0;
	double const ZERO  =  0.0;
}

// ----------------------------------------------------
class Model
{
public:

	Model( Logger * const logger );
	~Model(){};

	void Initialize(void);

	// assign values of these variables outside of class

	std::string   build;
	std::string   model_name;
	unsigned int  gcode;
	unsigned int  gene_min_length;
	unsigned int  order_cod;
	unsigned int  order_non;

	double  noncoding_duration_decay;
	double  coding_duration_decay;

	double  probability_N_in_noncoding;
	double  probability_N_in_coding;

	std::vector<double>  non;
	std::vector<double>  cod1;
	std::vector<double>  cod2;
	std::vector<double>  cod3;

	double  pATG;
	double  pGTG;
	double  pTTG;

	double  pTAA;
	double  pTAG;
	double  pTGA;

	double toModel;

	unsigned int ORF_start_marging;

	// assign values internaly - inside of class

	void ReserveSpace(void);
	
	// logodd:  log( cod[i]/non[i] )

	std::vector<double>  logodd_1;
	std::vector<double>  logodd_2;
	std::vector<double>  logodd_3;

	std::vector<double>  logodd_1_abs;
	std::vector<double>  logodd_2_abs;
	std::vector<double>  logodd_3_abs;
	
	// log probabiliti normalized on length:  log( p[i]/0.25^(order+1))

	// coding normalized

	std::vector<double>  logP_1_n;
	std::vector<double>  logP_2_n;
	std::vector<double>  logP_3_n;
	
	std::vector<double>  logP_1_abs_n;
	std::vector<double>  logP_2_abs_n;
	std::vector<double>  logP_3_abs_n;

	// noncoding normalized

	std::vector<double>  logP_N_n;
	std::vector<double>  logP_N_abs_n;

	// duration

	double GetDurationForORF(unsigned int const length);
	double GetIncompleteDurationForORF(unsigned int const length);

	// start 

	double logodd_ATG;
	double logodd_GTG;
	double logodd_TTG;

	// stop

	double logodd_TAA;
	double logodd_TAG;
	double logodd_TGA;

	// sites

	Site RBS;
	Site StartContentRBS;

	Site Promoter;
	Site StartContentPromoter;

	Site StartContent;
    
    // MGM sites
    Site RBS_A, RBS_B, RBS_C, RBS_D, RBS_X;
    Site SC_RBS_A, SC_RBS_B, SC_RBS_C, SC_RBS_D, SC_RBS_X;
    Site PROMOTER_C, PROMOTER_D;
    Site SC_PROMOTER_C, SC_PROMOTER_D;
    Site EUS;

	//

	std::string KL(void);

private:

	std::vector<double> ReduceOrderAbs( std::vector<double> const & arr , unsigned int const order_in, unsigned int const order_out );

	void VerifyProbabilityArray( std::vector<double> const & arr, unsigned int const order, std::string const & message );
	std::string IntToString( unsigned int index, const unsigned int order );
	void NormalizeArray( std::vector<double> & arr, std::string const & message );
	void CheckForNoZeroValues( std::vector<double> const & arr, unsigned int const order, std::string const & message );

	void CalculateAbsoluteCodVsNonRatio( std::vector<double> const & c, std::vector<double> const & n, std::vector<double> & target );
	void CalculateAbsoluteCodVsNonRatio( std::vector<double> const & c, std::vector<double> const & n, std::vector<double> const & nc, std::vector<double> & target );

	void AbsoluteToConditionalByLastPosition( std::vector<double> const & source, std::vector<double> & target );
	void DivideByConditionalNon( std::vector<double> const & nc, std::vector<double> & target );

	void SetLogOddForN( std::vector<double> & target, double const p );
	void Log(std::vector<double> & target);
	double LogRatio( double x, double y );

	void CalculateLogOdd(void);
	void CalculateStartStopLoggodd(void);

	void CalculateLogProb(void);

	unsigned int ReversCompIndex( unsigned int idx, unsigned int order );
	void ReverseCompNoncodingCounts( unsigned int order, std::vector<double> & arr );

	void InitializeSite( Site * ptr, bool with_dur, int norm_order, std::string const & label );

	std::string PrintLogOdds(void);
	std::string PrintLogOdds_abs(void);
	std::string PrintDuration(void);

	void InitiateLogoddDuration(unsigned int const reserved_length);
	void InitiateLogoddIncompleDuration(unsigned int const reserved_length);

	std::vector<double> logodd_duration;
	std::vector<double> logodd_incomplete_duration;

	double dur_c;
	double dur_r;

	Logger * const logger;
};
// ----------------------------------------------------
#endif // MODEL_2_H

