#ifndef SFM_DRIVER2_H__
#define SFM_DRIVER2_H__

#include "lm.h"
#include "dataf/camera.h"
#include "geometry_data.h"
#include <vector>
#include <matrix/vector.h>
#include "dataf/keypoint.h"
#include "dataf/dataf.h"

namespace bundle{
	
class GRTxyRefiner;

namespace{
// /*used to refine the global rotation*/
class GRRefiner{
private:
	double opts[LM_OPTS_SZ], info[LM_INFO_SZ];
	std::vector<std::vector<v3_t> > tracks;

private:
	static void modhorizon(double *p, double *x, int m, int n, void *data);
	static void modhorizonjac(double *p, double *jac, int m, int n, void *data);
	
private:
	double Run__(const GeometryData& data, double percent, double add_r, const double* tx, const double* ty, double* ggamma);
	
public:
	GRRefiner();
	~GRRefiner();
	
	double Run(const GeometryData& data, double percent, const double* tx, const double* ty, double* ggamma);
	
	friend class bundle::GRTxyRefiner;
};

/*used to refine the t_y*/
class TyRefiner{
private:
	double opts[LM_OPTS_SZ], info[LM_INFO_SZ];
	std::vector<std::vector<v3_t> > tracks;
	static int ankidx;

private:
	static void modhorizonty(double *p, double *x, int m, int n, void *data);
	
private:
	double Run__(const GeometryData& data, double percent, double gra, const double* tya, double* ty);		//gra: global rotation angle
	
public:
	TyRefiner();
	~TyRefiner();
	
	double Run(const GeometryData& data, double percent, double gra, double* ty);
	
	friend class bundle::GRTxyRefiner;
};

class XFormRefiner{
protected:
	static double xform(const std::vector<util::_point>& pts1, const std::vector<util::_point>& pts2, double xf[6]);	
};

class TxRefiner: public XFormRefiner{
private:
	static int ankidx;

private:
	void Run(const GeometryData& data, double percent, double add_r, const double* ty, double* tx);
	
public:
	
	friend class bundle::GRTxyRefiner;
};

class CentreRefiner: public XFormRefiner{
private:
	static int ankidx;
private:
	void Run(const GeometryData& data, double* tx, double* ty);
	
public:
	
	friend class bundle::GRTxyRefiner;
};

class SeriesRRefiner{
private:
	double opts[LM_OPTS_SZ], info[LM_INFO_SZ];
	std::vector<std::vector<v3_t> > tracks;

private:
	static void modhorizonty(double *p, double *x, int m, int n, void *data);
	
private:
	double Run__(const GeometryData& data, double percent, double gra, const double* tx, const double* ty, double* sr);		//gra: global rotation angle
	
public:
	SeriesRRefiner();
	~SeriesRRefiner();
	
	double Run(const bundle::GeometryData& data, double percent, double gra, const double* tx, const double* ty, double* sr);
	
	friend class bundle::GRTxyRefiner;
};

}

class GRTxyRefiner{
public:
	GRTxyRefiner();
	~GRTxyRefiner();
	
	void Run(const GeometryData& data, double percent, double* sgamma, double* tx, double* ty);
};

static void RecalculateXFormTxy(const double* sr, const double* tx, const double* ty, int n, double* ntx, double* nty);

static void RecalculateXFormTxy(double gr, const double* tx, const double* ty, int n, double* ntx, double* nty);

}

#endif
