#include "sfm_driver2.h"
#include <cmath>
#include <iostream>
#include "util/exception.h"
#include <cfloat>
#include "dataf/keypoint.h"
#include "geometry_data.h"
#include <unistd.h>
#include <numeric>
#include <lm.h>
#include <climits>
#include "matrix/matrix.h"

bundle::GRRefiner::~GRRefiner(){}

bundle::GRRefiner::GRRefiner()
{
	opts[0]=LM_INIT_MU; opts[1]=1E-15; opts[2]=1E-15; opts[3]=1E-20;opts[4]= LM_DIFF_DELTA;
}

void bundle::GRRefiner::modhorizon(double* p, double* x, int m, int n, void* data)
{
	const std::vector<std::vector<v3_t> >& track_data = *(std::vector<std::vector<v3_t> >*)data;
	
	double cosp = cos(p[0]);
	double sinp = sin(p[0]);
	
	for(int i = 0; i < track_data.size(); i++){
		const std::vector<v3_t>& track = track_data[i];
		std::vector<float> new_y;
		for(int j = 0; j < track.size(); j++){
			const v3_t& key = track[j];
			new_y.push_back(-sinp*key.p[0]+cosp*key.p[1]);
		}
		float avg_ny = std::accumulate(new_y.begin(), new_y.end(), 0.0f)/new_y.size();
		
		x[i] = 0;
		for(int j = 0; j < new_y.size(); j++){
			x[i] += fabs(new_y[j]-avg_ny);
		}
	}
}

void bundle::GRRefiner::modhorizonjac(double* p, double* jac, int m, int n, void* data)
{
	const std::vector<std::vector<v3_t> >& track_data = *(std::vector<std::vector<v3_t> >*)data;
	
	double cosp = cos(p[0]);
	double sinp = sin(p[0]);
	
	for(int i = 0; i < track_data.size(); i++){
		const std::vector<v3_t>& track = track_data[i];
		std::vector<float> jac_y;
		for(int j = 0; j < track.size(); j++){
			const v3_t& key = track[j];
			jac_y.push_back(-cosp*key.p[0]-sinp*key.p[1]);
		}
		float avg_ny = std::accumulate(jac_y.begin(), jac_y.end(), 0.0f)/jac_y.size();
		
		jac[i] = 0;
		for(int j = 0; j < jac_y.size(); j++){
			jac[i] += fabs(jac_y[j]-avg_ny);
		}
	}
}

double bundle::GRRefiner::Run__(const bundle::GeometryData& data, double percent, double add_r, const double* tx, const double* ty, double* ggamma)
{	
	const std::vector<TrackData>& track_data = data.m_track_data;
	tracks.clear();
	
	int length_thre = percent*data.NumImages();
	double cosar = cos(add_r);
	double sinar = sin(add_r);
	
	for(int i = 0; i < track_data.size(); i++){
		const TrackData& track = track_data[i];
		
		if(track.m_views.size() < length_thre){
			continue;
		}
		
		std::vector<v3_t> ntrc;
	
		for(int j = 0; j < track.m_views.size(); j++){
			const Keypoint& key = data.GetKey(track.m_views[j].first, track.m_views[j].second);
			ntrc.push_back(v3_new(cosar*key.m_x+sinar*key.m_y , -sinar*key.m_x+cosar*key.m_y, track.m_views[j].first));
		}
		
		tracks.push_back(ntrc);
	}
	
	if(tx){
		for(int i = 0; i < tracks.size(); i++){
			std::vector<v3_t>& track = tracks[i];
			for(int j = 0; j < track.size(); j++){
// 				std::cout<<track[j].p[0]+ tx[(int)track[j].p[2]]<<" "<<track[j].p[1]+ ty[(int)track[j].p[2]]<<std::endl;
				track[j].p[0] += tx[(int)track[j].p[2]];
			}
// 			std::cout<<std::endl;
		}
	}
	
	if(ty){
		for(int i = 0; i < tracks.size(); i++){
			std::vector<v3_t>& track = tracks[i];
			for(int j = 0; j < track.size(); j++){
				track[j].p[1] += ty[(int)track[j].p[2]];
			}
		}
	}
	
	int m = 1, n = tracks.size();
	double* x = new double[n];
	
	memset(x, 0, sizeof(double)*n);
	int ret = dlevmar_der(modhorizon, modhorizonjac, ggamma, x, m, n, 1000, opts, info, NULL, NULL, &tracks);
	
	modhorizon(ggamma, x, m, n, &tracks);
	
	double res = std::accumulate(x, x+n, 0.0f)/n;
	
	delete [] x;
	
	return res;
}

static double FormatAngle(double angle)
{
	while(angle > 1.5707963267949 || angle < -1.5707963267949){
		if(angle > 1.5707963267949){
			angle -= M_PI;
		}
		else if(angle < -1.5707963267949){
			angle += M_PI;
		}
	}
	
	return angle;
}

double bundle::GRRefiner::Run(const bundle::GeometryData& data, double percent, const double* tx, const double* ty, double* ggamma)
{
	float best_angle;
	float min_res = INT_MAX;
	
	for(float r = -1.5707963267949; r <= 1.5707963267949; r += 0.0872664625997165){
		
		double cr = cos(r), sr = sin(r);
		double ntx[data.NumImages()], nty[data.NumImages()];
		
		for(int i = 0; i < data.NumImages(); i++){
			ntx[i] = cr*tx[i]+sr*ty[i];
			nty[i] = -sr*tx[i]+cr*ty[i];
		}
		
		*ggamma = 0;
		double res = Run__(data, percent, r, ntx, nty, ggamma);
		
		*ggamma = FormatAngle(*ggamma+r);
		
		if(res < min_res){
			min_res = res;
			best_angle = *ggamma;
		}
		
// 		std::cout<<r/M_PI*180<<"\t"<<*ggamma/M_PI*180<<"\t"<<res<<std::endl;
	}
	
	*ggamma = best_angle;
	return min_res;
}

bundle::SeriesRRefiner::~SeriesRRefiner(){}

bundle::SeriesRRefiner::SeriesRRefiner()
{
	opts[0]=LM_INIT_MU; opts[1]=1E-15; opts[2]=1E-15; opts[3]=1E-20;opts[4]= LM_DIFF_DELTA;
}

void bundle::SeriesRRefiner::modhorizonty(double* p, double* x, int m, int n, void* data)
{
	const std::vector<std::vector<v3_t> >& track_data = *(std::vector<std::vector<v3_t> >*)data;
	
	double cosv[m];
	double sinv[m];
	
	for(int i = 0; i < m; i++){
		cosv[i] = cos(p[i]);
		sinv[i] = sin(p[i]);
	}
	
	memset(x, 0, sizeof(double)*n);
	
	for(int i = 0; i < track_data.size(); i++){
		const std::vector<v3_t>& track = track_data[i];
		std::vector<float> new_y;
		for(int j = 0; j < track.size(); j++){
			const v3_t& key = track[j];
			int s_idx = (int)key.p[2];
			new_y.push_back(-sinv[s_idx]*key.p[0]+cosv[s_idx]*key.p[1]);
		}
		float avg_ny = std::accumulate(new_y.begin(), new_y.end(), 0.0f)/new_y.size();
		
		for(int j = 0; j < track.size(); j++){
			x[(int)(track[j].p[2])] += fabs(new_y[j]-avg_ny);
		}
	}
}

double bundle::SeriesRRefiner::Run__(const bundle::GeometryData& data, double percent, double gra, const double* tx, const double* ty, double* sr)
{
	const std::vector<TrackData>& track_data = data.m_track_data;
	tracks.clear();
	
	int length_thre = percent*data.NumImages();
	double cosar = cos(gra);
	double sinar = sin(gra);
	
	for(int i = 0; i < track_data.size(); i++){
		const TrackData& track = track_data[i];
		
		if(track.m_views.size() < length_thre){
			continue;
		}
		
		std::vector<v3_t> ntrc;
	
		for(int j = 0; j < track.m_views.size(); j++){
			const Keypoint& key = data.GetKey(track.m_views[j].first, track.m_views[j].second);
			ntrc.push_back(v3_new(cosar*key.m_x+sinar*key.m_y , -sinar*key.m_x+cosar*key.m_y, track.m_views[j].first));
		}
		
		tracks.push_back(ntrc);
	}
	
	if(tx){
		for(int i = 0; i < tracks.size(); i++){
			std::vector<v3_t>& track = tracks[i];
			for(int j = 0; j < track.size(); j++){
				track[j].p[0] += tx[(int)track[j].p[2]];
			}
		}
	}
	
	if(ty){
		for(int i = 0; i < tracks.size(); i++){
			std::vector<v3_t>& track = tracks[i];
			for(int j = 0; j < track.size(); j++){
				track[j].p[1] += ty[(int)track[j].p[2]];
			}
		}
	}
	
	int m = data.NumImages(), n = data.NumImages();
	double* x = new double[n];
	memset(x, 0, sizeof(double)*n);
	
	int ret = dlevmar_dif(modhorizonty, sr, x, m, n, 1000, opts, info, NULL, NULL, &tracks);
	
	modhorizonty(sr, x, m, n, &tracks);
	
	double res = std::accumulate(x, x+n, 0.0f)/n;
	
	delete [] x;
	
	return res;
}

double bundle::SeriesRRefiner::Run(const bundle::GeometryData& data, double percent, double gra, const double* tx, const double* ty, double* sr)
{
	return Run__(data, percent, gra, tx, ty, sr);
}

bundle::TyRefiner::~TyRefiner(){}

bundle::TyRefiner::TyRefiner()
{
	opts[0]=LM_INIT_MU; opts[1]=1E-15; opts[2]=1E-15; opts[3]=1E-20;opts[4]= LM_DIFF_DELTA;
}

int bundle::TyRefiner::ankidx = 0;

void bundle::TyRefiner::modhorizonty(double* p, double* x, int m, int n, void* data)
{
	const std::vector<std::vector<v3_t> >& track_data = *(std::vector<std::vector<v3_t> >*)data;
	p[ankidx] = 0;
	memset(x, 0, sizeof(double)*n);
	
	for(int i = 0; i < track_data.size(); i++){
		const std::vector<v3_t>& track = track_data[i];
		std::vector<float> new_y;
		for(int j = 0; j < track.size(); j++){
			const v3_t& key = track[j];
			new_y.push_back(key.p[1]+p[(int)key.p[2]]);
		}
		float avg_ny = std::accumulate(new_y.begin(), new_y.end(), 0.0f)/new_y.size();
		
		for(int j = 0; j < track.size(); j++){
			x[(int)(track[j].p[2])] += fabs(new_y[j]-avg_ny);
		}
	}
}

double bundle::TyRefiner::Run__(const bundle::GeometryData& data, double percent, double gra, const double* tya, double* ty)
{
	const std::vector<TrackData>& track_data = data.m_track_data;
	tracks.clear();
	
	int length_thre = percent*data.NumImages();
	double cosar = cos(gra);
	double sinar = sin(gra);
	
	float min_rot = 360;
	for(int i = 0; i < data.NumImages(); i++){
		float absang = fabs(data.OriAngles(i));
		if(absang < min_rot){
			min_rot = absang;
			ankidx = i;
		}
	}
	
	for(int i = 0; i < track_data.size(); i++){
		const TrackData& track = track_data[i];
		
		if(track.m_views.size() < length_thre){
			continue;
		}
		
		std::vector<v3_t> ntrc;
	
		for(int j = 0; j < track.m_views.size(); j++){
			const Keypoint& key = data.GetKey(track.m_views[j].first, track.m_views[j].second);
			ntrc.push_back(v3_new(cosar*key.m_x+sinar*key.m_y , -sinar*key.m_x+cosar*key.m_y, track.m_views[j].first));
		}
		
		tracks.push_back(ntrc);
	}
	
	double cptya[data.NumImages()];
	
	if(tya){
		for(int i = 0; i < tracks.size(); i++){
			std::vector<v3_t>& track = tracks[i];
			for(int j = 0; j < track.size(); j++){
				track[j].p[1] += tya[(int)track[j].p[2]];
			}
		}
		memcpy(cptya, tya, sizeof(double)*data.NumImages());
	}
	
	int m = data.NumImages(), n = data.NumImages();
	double* x = new double[n];
	memset(x, 0, sizeof(double)*n);
	
	int ret = dlevmar_dif(modhorizonty, ty, x, m, n, 1000, opts, info, NULL, NULL, &tracks);
	
	modhorizonty(ty, x, m, n, &tracks);
	
	double res = std::accumulate(x, x+n, 0.0f)/n;
	
	delete [] x;
	
	if(tya){
		for(int i = 0; i < data.NumImages(); i++){
			ty[i] += cptya[i];
		}
	}
	
	return res;
}

double bundle::TyRefiner::Run(const bundle::GeometryData& data, double percent, double gra, double* ty)
{
	return Run__(data, percent, gra, NULL, ty);
}

int bundle::TxRefiner::ankidx = 0;

/**xf: a, b, c; d, e, f**/
double bundle::XFormRefiner::xform(const std::vector<util::_point>& pts1, const std::vector<util::_point>& pts2, double xf[6])
{
	assert(pts1.size() == pts2.size());
	int num_eqs = 2*pts1.size();
    int num_vars = 6;

    double *As = new double[num_eqs*num_vars];
    double *bs = new double[num_eqs];

    for(int i = 0; i < pts1.size(); i++){
        double* A = As+12*i;
        double* b = bs+2*i;
        A[0] = pts1[i].x; A[1] = pts1[i].y; A[2] = 1; A[3] = 0; A[4] = 0; A[5] = 0;
		A[6] = 0; A[7] = 0; A[8] = 0; A[9] = pts1[i].x; A[10] = pts1[i].y; A[11] = 1;
		
        
        b[0] = pts2[i].x;
        b[1] = pts2[i].y;
    }
    
    dgelsy_driver(As, bs, xf, num_eqs, num_vars, 1);
	
	double res = 0;
	for(int i = 0; i < pts1.size(); i++){
		float deltax = pts1[i].x*xf[0]+pts1[i].y*xf[1]+xf[2];
		float deltay = pts1[i].x*xf[3]+pts1[i].y*xf[4]+xf[5];
        
		deltax -= pts2[i].x;
		deltay -= pts2[i].y;
		
		res += deltax*deltax+deltay*deltay;
    }
    
    res = sqrt(res)/pts1.size();
    
    delete [] As;
    delete [] bs;

    return res;
}

void bundle::TxRefiner::Run(const bundle::GeometryData& data, double percent, double add_r, const double* ty, double* tx)
{	
	double cosar = cos(add_r);
	double sinar = sin(add_r);
	
	float min_rot = 360;
	for(int i = 0; i < data.NumImages(); i++){
		float absang = fabs(data.OriAngles(i));
		if(absang < min_rot){
			min_rot = absang;
			ankidx = i;
		}
	}
	
	const MatchTable& matchtable = data.m_matches;
	
	for(int i = 0; i < data.NumImages()-1; i++){
		const std::vector<KeypointMatch>& matchidx = matchtable.GetMatchList(MatchIndex(i, i+1));
		std::vector<util::_point> points1;
		std::vector<util::_point> points2;
		
		for(int j = 0; j < matchidx.size(); j++){
			Keypoint key; util::_point pt;
			
			key = data.GetKey(i, matchidx[j].m_idx1);
			pt.x = cosar*key.m_x+sinar*key.m_y; 
			pt.y = -sinar*key.m_x+cosar*key.m_y+ty[i];
			points1.push_back(pt);
			
			key = data.GetKey(i+1, matchidx[j].m_idx2);
			pt.x = cosar*key.m_x+sinar*key.m_y; 
			pt.y = -sinar*key.m_x+cosar*key.m_y+ty[i+1];
			points2.push_back(pt);
		}
		
		double xf[6];
		double res = xform(points1, points2, xf);
// 		util::_point centre;
// 		centre.y = ((xf[0]-1)*xf[5]-xf[2]*xf[3])/(xf[3]*xf[1]-(xf[0]-1)*(xf[4]-1));
// 		centre.x = (xf[2]+xf[1]*centre.y)/(1-xf[0]);
		
		tx[i] = xf[2];
	}
	
	for(int i = ankidx-2; i >= 0; i--){
		tx[i] += tx[i+1];
 	}
	
	for(int i = data.NumImages()-1; i >=ankidx+1; i--){
		tx[i] = -tx[i-1];
	}
	
	for(int i = ankidx+2; i < data.NumImages(); i++){
		tx[i] += tx[i-1];
	}
	
	tx[ankidx] = 0;
}

int bundle::CentreRefiner::ankidx = 0;

void bundle::CentreRefiner::Run(const bundle::GeometryData& data, double* tx, double* ty)
{
	float min_rot = 360;
	for(int i = 0; i < data.NumImages(); i++){
		float absang = fabs(data.OriAngles(i));
		if(absang < min_rot){
			min_rot = absang;
			ankidx = i;
		}
	}
	
	const MatchTable& matchtable = data.m_matches;
	
	double xfv[data.NumImages()-1][6];
	std::vector<util::_point> centre(data.NumImages());
	
	for(int i = 0; i < data.NumImages()-1; i++){
		const std::vector<KeypointMatch>& matchidx = matchtable.GetMatchList(MatchIndex(i, i+1));
		std::vector<util::_point> points1;
		std::vector<util::_point> points2;
		
		for(int j = 0; j < matchidx.size(); j++){
			Keypoint key; util::_point pt;
			
			key = data.GetKey(i, matchidx[j].m_idx1);
			pt.x = key.m_x; 
			pt.y = key.m_y;
			points1.push_back(pt);
			
			key = data.GetKey(i+1, matchidx[j].m_idx2);
			pt.x = key.m_x; 
			pt.y = key.m_y;
			points2.push_back(pt);
		}
		
		double res = xform(points1, points2, xfv[i]);
// 		std::cout<<"xform res: "<<res<<std::endl;
	}
	
	centre[ankidx].x = 0; centre[ankidx].y = 0;
	
	for(int i = ankidx+1; i < data.NumImages(); i++){
		double* xf = xfv[i-1];
		
		centre[i].x = xf[0]*centre[i-1].x+xf[1]*centre[i-1].y+xf[2];
		centre[i].y = xf[3]*centre[i-1].x+xf[4]*centre[i-1].y+xf[5];
	}
	
	for(int i = ankidx-1; i >=0; i--){
		double* xf = xfv[i];
		double A[] = {xf[0], xf[1], xf[3], xf[4]};
		double A_inv[4];
		matrix_invert(2, A, A_inv);
		
		double nx = centre[i+1].x-xf[2], ny = centre[i+1].y-xf[5];
		
		centre[i].x = A_inv[0]*nx+A_inv[1]*ny;
		centre[i].y = A_inv[2]*nx+A_inv[3]*ny;
	}
	
	for(int i = 0; i < data.NumImages(); i++){
		tx[i] = -centre[i].x;
		ty[i] = -centre[i].y;
		
// 		std::cout<<"xy: "<<tx[i]<<"\t"<<ty[i]<<std::endl;
	}
}

bundle::GRTxyRefiner::GRTxyRefiner(){}

bundle::GRTxyRefiner::~GRTxyRefiner(){}

void bundle::GRTxyRefiner::Run(const bundle::GeometryData& data, double percent, double* sgamma, double* tx, double* ty)
{
	
	GRRefiner grrefiner;
	TyRefiner tyrefiner;
	TxRefiner txrefiner;
	CentreRefiner crefiner;
	SeriesRRefiner srrefiner;
	double res;
	
	double ggamma[1];
	
	crefiner.Run(data, tx, ty);
	
	res = grrefiner.Run(data, percent, tx, ty, ggamma);
	RecalculateXFormTxy(*ggamma, tx, ty, data.NumImages(), tx, ty);
	res = tyrefiner.Run(data, percent, *ggamma, ty);
	
	std::cout<<*ggamma/M_PI*180<<"\t"<<res<<std::endl;
	
	double res_old = DBL_MAX;
	int count = 0;
	double gmres = DBL_MAX, subgmres = DBL_MAX;
	
	while(true){
		double gammp;  //, tymp[data.NumImages()], txmp[data.NumImages()];
		gammp = *ggamma; *ggamma = 0;
// 		memcpy(tymp, ty, sizeof(double)*data.NumImages());
// 		memcpy(txmp, tx, sizeof(double)*data.NumImages());
		
		float best_angle;
// 		double added_angle;
		float min_res = INT_MAX;
		for(float r = -0.0872664625997165; r <= 0.0872664625997165; r += 0.0174532925199433){
			*ggamma = r;
			res = grrefiner.Run__(data, percent, gammp, tx, ty, ggamma);
			
			if(res < min_res){
				min_res = res;
				*ggamma = FormatAngle(*ggamma+gammp);
				best_angle = *ggamma;
			}
		}
		
		*ggamma = best_angle;
		
		for(int i = 0; i < data.NumImages(); i++){
			ty[i] += ((rand()%100)-50.f)/500;
		}
		
		res = tyrefiner.Run(data, percent, *ggamma, ty);
		
		if(res < gmres){
			subgmres = gmres;
			gmres = res;
		}
		else if(res < subgmres){
			subgmres = res;
		}
	
		std::cout<<*ggamma/M_PI*180<<"\t"<<res<<std::endl;
		
// 		for(int i = 0; i < data.NumImages(); i++){
// 			std::cout<<ty[i]<<" ";
// 		}
// 		std::cout<<std::endl;
		
		if(subgmres-gmres < 0.3 && res-gmres < 0.3 && fabs(res-gmres) < 0.3 && fabs(subgmres-gmres) < 0.3){
			break;
		}
		
		if(fabs(res-res_old) < 0.3){
			if(count++ > 5){
				break;
			}
		}
		else{
			count = 0;
		}
		
		res_old = res;
		
		txrefiner.Run(data, 0, *ggamma, ty, tx);
	}

	txrefiner.Run(data, 0, *ggamma, ty, tx);
	srrefiner.Run(data, percent, *ggamma, tx, ty, sgamma);
	
	for(int i = 0; i < data.NumImages(); i++){
		sgamma[i] += *ggamma;
// 		std::cout<<sgamma[i]/M_PI*180<<std::endl;
	}
// 	for(int i = 0; i < data.NumImages(); i++){
// 		std::cout<<tx[i]<<"\t"<<ty[i]<<std::endl;
// 	}
	
	return;
}

void bundle::RecalculateXFormTxy(double gr, const double* tx, const double* ty, int n, double* ntx, double* nty)
{
	double sr[n];
	for(int i = 0; i < n; i++){
		sr[i] = gr;
	}
	
	RecalculateXFormTxy(sr, tx, ty, n, ntx, nty);
}

void bundle::RecalculateXFormTxy(const double* sr, const double* tx, const double* ty, int n, double* ntx, double* nty)
{
	for(int i = 0; i < n; i++){
		double cosr = cos(sr[i]), sinr = sin(sr[i]);
		double x = tx[i], y = ty[i];
		
		ntx[i] = cosr*x+sinr*y;
		nty[i] = -sinr*x+cosr*y;
	}
}
