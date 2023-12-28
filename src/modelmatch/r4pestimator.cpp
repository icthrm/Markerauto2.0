#include "r4pestimator.h"
#include "searchtree.h"
// #include "util/matrix.h"
// #include "matrix/matrix.h"

// Ran4PEstimator::Ran4PEstimator(std::vector< util::point2d >& __fids1, std::vector< util::point2d >& __fids2):fids1(__fids1), fids2(__fids2)
// {
// 	srand((unsigned)time( NULL ));
// 	cos_beta1 = INT_MAX; cos_beta2 = INT_MAX;
// }

Ran4PEstimator::Ran4PEstimator(std::vector< util::point2d >& __fids1, std::vector< util::point2d >& __fids2, bool method):fids1(__fids1), fids2(__fids2)
{
	if(method)
	{
		nPts = -1;
		nk = KDTREE_NK;
		eps = KDTREE_EPS;
		quad_svd_safe_tol = QUAD_SVD_SAFE_TOL;
		graph_baseline_safe_tol = GRAPH_BASE_SAFE_TOL_old;
		queryPt = NULL;
		dataPts = NULL;
		nnIdx = NULL;
		dists = NULL;
		kdTree = NULL;
		srand((unsigned)time( NULL ));
	}
	else{
		srand((unsigned)time( NULL ));
		cos_beta1 = INT_MAX; cos_beta2 = INT_MAX;
	}
}

Ran4PEstimator::~Ran4PEstimator(){}

void Ran4PEstimator::PointPair::FractionPoint(float fraction, util::point2d& fract)
{
	fract.x = p1->x+fraction*p1top2.x;
	fract.y = p1->y+fraction*p1top2.y;
}

#define PI_180				0.0174532925199433

void Ran4PEstimator::SetSupplementData(double beta1, double beta2)
{
	cos_beta1 = cos(beta1*PI_180); cos_beta2 = cos(beta2*PI_180);
}

void Ran4PEstimator::PairDistanceMeanStd(const std::vector<PointPair>& pointpairs, float* dist_mean, float* dist_std)
{
	*dist_mean = 0; 
	*dist_std = 0;
	
	for(long i = pointpairs.size(); i--; ){
		*dist_std += pointpairs[i].squared_length;
		*dist_mean += sqrt(pointpairs[i].squared_length);
	}
	
	*dist_std /= pointpairs.size();
	*dist_mean /= pointpairs.size();
	*dist_std -= *dist_mean* *dist_mean;
	*dist_std = sqrt(*dist_std);
}

void Ran4PEstimator::GeneratePointPairs(const std::vector< util::point2d >& fids, std::vector< Ran4PEstimator::PointPair >& pointpairs)
{
	for(int i = 0; i < fids.size(); i++){
		for(int j = i+1; j < fids.size(); j++){
			PointPair ppair(i, j, fids[i], fids[j]);
			pointpairs.push_back(ppair);
		}
	}
}

void Ran4PEstimator::FilterPointPairs(float max_dist, float min_dist, std::vector< Ran4PEstimator::PointPair >& pointpairs)
{
	float max_dist2 = max_dist*max_dist;
	float min_dist2 = min_dist*min_dist;
	
	std::vector< Ran4PEstimator::PointPair > filtered;
// 	std::cout<<pointpairs.size()<<std::endl;
	
	for(int i = pointpairs.size(); i--;){
		if(pointpairs[i].squared_length < min_dist2 || pointpairs[i].squared_length > max_dist2){
			continue;
		}
		filtered.push_back(pointpairs[i]);
	}
	pointpairs = filtered;
	
// 	std::cout<<pointpairs.size()<<std::endl;
}

static void GenerateRandomIndex(int size, int* idxs, int idx_size){
	idxs[0] = rand()%size;
	int idx = 1;
	while(true){
generatebeginrqervfe124334321:
		int tmp = rand()%size;
		for(int i = 0; i < idx; i++){
			if(tmp == idxs[i]){
				goto generatebeginrqervfe124334321;
			}
		}
		idxs[idx] = tmp;
		if(++idx >= idx_size){
			break;
		}
	}
}

static inline bool IfTwoLineIntersection(const util::point2d& p1, const util::point2d& p2, const util::point2d& q1, const util::point2d& q2)
{
	if(!(min(p1.x,p2.x) <= max(q1.x,q2.x) && min(q1.x,q2.x) <= max(p1.x,p2.x) && min(p1.y,p2.y) <= max(q1.y,q2.y) && min(q1.y,q2.y) <= max(p1.y,p2.y)))
		return false;
	
	if(!(((q1.x-p1.x)*(q1.y-q2.y)-(q1.y-p1.y)*(q1.x-q2.x))*((q1.x-p2.x)*(q1.y-q2.y)-(q1.y-p2.y)*(q1.x-q2.x)) < 0 &&
		((p1.x-q1.x)*(p1.y-p2.y)-(p1.y-q1.y)*(p1.x-p2.x))*((p1.x-q2.x)*(p1.y-p2.y)-(p1.y-q2.y)*( p1.x-p2.x)) < 0))
		return false;
	
	return true;
}

void Ran4PEstimator::SetPointSetDistErrorTOL(double psdetol)
{
	pts_dist_tol = psdetol;
}

int Ran4PEstimator::NK() const
{
	return nk;
}

double Ran4PEstimator::Eps() const
{
	return eps;
}


static inline void GetIntersectionP(const util::point2d& p1, const util::point2d& p2, const util::point2d& p3, const util::point2d& p4, util::point2d& c)
{
	int m = 2, n = 2, nrhs = 1;
	float As[4], *bs = c.v;
	As[0] = p2.y-p1.y; As[1] = p4.y-p3.y;
	As[2] = -p2.x+p1.x; As[3] = -p4.x+p3.x;
	
	bs[0] = (p2.y-p1.y)*p1.x-(p2.x-p1.x)*p1.y;
	bs[1] = (p4.y-p3.y)*p3.x-(p4.x-p3.x)*p3.y;
	
	int lda = m, ldb = m;
	int jpvt[n]; memset(jpvt, 0, sizeof(int)*n);
	float rcond = -1.0;
	int rank; /* Output */
	int lwork = -1;
	float test;//, *work; //= (float*)malloc(sizeof(float) * 1);
	int info;

	/* Query to find a good size for the work array */
	sgelsy_(&m, &n, &nrhs, As, &lda, bs, &ldb, jpvt, &rcond, &rank, &test, &lwork, &info);

	lwork = (int) test;//work[0];
	/* printf("Work size: %d\n", lwork); */
// 	free(work);
	float work[lwork];// = (float*)malloc(sizeof(float) * lwork);

	/* Make the FORTRAN call */
	sgelsy_(&m, &n, &nrhs, As, &lda, bs, &ldb, jpvt, &rcond, &rank, work, &lwork, &info);
}

void Ran4PEstimator::PickRandomQuadrilateral(const std::vector< Ran4PEstimator::PointPair >& ppairs, fpcs::Quadrilateral& quad, float* invariant1, float* invariant2)
{
	bool picked = false;
	while(!picked){
		int idx;
		GenerateRandomIndex(ppairs.size(), &idx, 1);
		
		for(int i = ppairs.size(); i--;){
			if(i == idx){
				continue;
			}
			
			if(!IfTwoLineIntersection(*(ppairs[idx].p1), *(ppairs[idx].p2), *(ppairs[i].p1), *(ppairs[i].p2))){
				continue;
			}
			
			util::point2d c;
			GetIntersectionP(*(ppairs[idx].p1), *(ppairs[idx].p2), *(ppairs[i].p1), *(ppairs[i].p2), c);
			
			*invariant1 = sqrt(util::L2(*(ppairs[idx].p1), c)/ppairs[idx].squared_length);
			*invariant2 = sqrt(util::L2(*(ppairs[i].p1), c)/ppairs[i].squared_length);
			
			quad[0] = ppairs[idx].idx1;
			quad[1] = ppairs[idx].idx2;
			quad[2] = ppairs[i].idx1;
			quad[3] = ppairs[i].idx2;
			
			quad.S = util::CalTriangleArea(ppairs[idx].P1(), ppairs[idx].P2(), ppairs[i].P1())+util::CalTriangleArea(ppairs[idx].P1(), ppairs[idx].P2(), ppairs[i].P2());
			
			picked = true;
			break;
		}
	}
}

void Ran4PEstimator::SearchQuadrilateralPairs(const fpcs::Quadrilateral& quad, float invariant1, float invariant2, float area_diff_tol, float dist_err_tol, std::vector< QuadPair >& pairs)
{
	pairs.clear();
	
	std::vector<util::point2d> e11(ppairs2.size());
	std::vector<util::point2d> e12(ppairs2.size());
	
	std::vector<util::point2d> e21(ppairs2.size());
	std::vector<util::point2d> e22(ppairs2.size());
	
// 	std::cout<<ppairs2.size()<<std::endl;
	
	util::point2d fract;
	float invariant12 = 1 - invariant1;
	float invariant22 = 1 - invariant2;
	for(int i = ppairs2.size(); i--; ){
		ppairs2[i].FractionPoint(invariant1, fract);
		e11[i] = (fract);
		ppairs2[i].FractionPoint(invariant12, fract);
		e12[i] = (fract);
		
		ppairs2[i].FractionPoint(invariant2, fract);
		e21[i] = (fract);
		ppairs2[i].FractionPoint(invariant22, fract);
		e22[i] = (fract);
	}
	
	SearchTree e21tree(e21, false);
	SearchTree e22tree(e22, false);
	
	float dist_err_thre = dist_err_tol*dist_err_tol;
	
	int idx; double distance;
	int* idices = &idx;
	double* distances = &distance;
	
	double area_ratio_l_boundary = (1-area_diff_tol)*quad.S*cos_beta2;
	double area_ratio_r_boundary = (1+area_diff_tol)*quad.S*cos_beta2;
	for(int i = e11.size(); i--;){
		e21tree.GetKNearestPoints(e11[i], 1, idices, distances);
		if(distances[0] < dist_err_thre){
			double S2 = util::CalTriangleArea(ppairs2[i].P1(), ppairs2[i].P2(), ppairs2[idices[0]].P1())+util::CalTriangleArea(ppairs2[i].P1(), ppairs2[i].P2(), ppairs2[idices[0]].P2());
			S2 = S2*cos_beta1;
			if(S2 < area_ratio_l_boundary || S2 > area_ratio_r_boundary){
				continue;
			}
			
			if(i == idices[0]){
				continue;
			}
			
			QuadPair qpair;
			qpair.first = quad;
			qpair.second.v[0] = ppairs2[i].idx1;
			qpair.second.v[1] = ppairs2[i].idx2;
			qpair.second.v[2] = ppairs2[idices[0]].idx1;
			qpair.second.v[3] = ppairs2[idices[0]].idx2;
			pairs.push_back(qpair);
		}
	}
	
	for(int i = e11.size(); i--;){
		e22tree.GetKNearestPoints(e11[i], 1, idices, distances);
		if(distances[0] < dist_err_thre){
			double S2 = util::CalTriangleArea(ppairs2[i].P1(), ppairs2[i].P2(), ppairs2[idices[0]].P1())+util::CalTriangleArea(ppairs2[i].P1(), ppairs2[i].P2(), ppairs2[idices[0]].P2());
			S2 = S2*cos_beta1;
			if(S2 < area_ratio_l_boundary || S2 > area_ratio_r_boundary){
				continue;
			}
			
			if(i == idices[0]){
				continue;
			}
			
			QuadPair qpair;
			qpair.first = quad;
			qpair.second.v[0] = ppairs2[i].idx1;
			qpair.second.v[1] = ppairs2[i].idx2;
			qpair.second.v[2] = ppairs2[idices[0]].idx2;
			qpair.second.v[3] = ppairs2[idices[0]].idx1;
// 			std::cout<<"S:"<<qpair.second.v[0]<<std::endl;
			pairs.push_back(qpair);
		}
	}
// 					std::cout<<"S:"<<pairs.size()<<std::endl;
	for(int i = e12.size(); i--;){
		e21tree.GetKNearestPoints(e12[i], 1, idices, distances);
		if(distances[0] < dist_err_thre){
			double S2 = util::CalTriangleArea(ppairs2[i].P1(), ppairs2[i].P2(), ppairs2[idices[0]].P1())+util::CalTriangleArea(ppairs2[i].P1(), ppairs2[i].P2(), ppairs2[idices[0]].P2());
			S2 = S2*cos_beta1;
			if(S2 < area_ratio_l_boundary || S2 > area_ratio_r_boundary){
				continue;
			}
			
			if(i == idices[0]){
				continue;
			}
			
			QuadPair qpair;
			qpair.first = quad;
			qpair.second.v[0] = ppairs2[i].idx2;
			qpair.second.v[1] = ppairs2[i].idx1;
			qpair.second.v[2] = ppairs2[idices[0]].idx1;
			qpair.second.v[3] = ppairs2[idices[0]].idx2;
			pairs.push_back(qpair);
		}
	}
	
	for(int i = e12.size(); i--;){
		e22tree.GetKNearestPoints(e12[i], 1, idices, distances);
		if(distances[0] < dist_err_thre){
			double S2 = util::CalTriangleArea(ppairs2[i].P1(), ppairs2[i].P2(), ppairs2[idices[0]].P1())+util::CalTriangleArea(ppairs2[i].P1(), ppairs2[i].P2(), ppairs2[idices[0]].P2());
			S2 = S2*cos_beta1;
			if(S2 < area_ratio_l_boundary || S2 > area_ratio_r_boundary){
				continue;
			}
			
			if(i == idices[0]){
				continue;
			}
			
			QuadPair qpair;
			qpair.first = quad;
			qpair.second.v[0] = ppairs2[i].idx2;
			qpair.second.v[1] = ppairs2[i].idx1;
			qpair.second.v[2] = ppairs2[idices[0]].idx2;
			qpair.second.v[3] = ppairs2[idices[0]].idx1;
			pairs.push_back(qpair);
		}
	}
}

Ran4PEstimator* Ran4PEstimator::HIterator::r4per = NULL;
cv::Mat Ran4PEstimator::HIterator::HH;
std::vector<std::pair<int, int> >* Ran4PEstimator::HIterator::matches = NULL;

//H = (a, b, c; d, e, f; g, h, i); x = {a, b, c, d, e, f, g, h, i}
void Ran4PEstimator::HIterator::lmderfun(double* p, double* x, int m, int n, void* data)
{
	std::vector<util::point2d>* pts = (std::vector<util::point2d>*)data;
	int count = 0;
	for(int i = 0; i < n; i += 2){
		float px = (*pts)[count].x, py = (*pts)[count].y;
		count++;
		x[i] = p[0]*px+p[1]*py+p[4];
		x[i+1] = p[2]*px+p[3]*py+p[5];
	}
}

void Ran4PEstimator::HIterator::lmderjac(double* p, double* jac, int m, int n, void* data)
{
	std::vector<util::point2d>* pts = (std::vector<util::point2d>*)data;
	int count = 0;
	int j = 0;
	for(int i = 0; i < n; i += 2){
		float px = (*pts)[count].x, py = (*pts)[count].y;
		count++;
		
		jac[j++] = px; jac[j++] = py; jac[j++] = 0; jac[j++] = 0; jac[j++] = 1; jac[j++] = 0;
		jac[j++] = 0; jac[j++] = 0; jac[j++] = px; jac[j++] = py; jac[j++] = 0; jac[j++] = 1;
	}	
}

void Ran4PEstimator::HIterator::SetInitialInput(Ran4PEstimator* __r4per, cv::Mat& __HH, std::vector< std::pair< int, int > >* __matches)
{
	r4per = __r4per;
	HH = __HH;
	matches = __matches;
}

bool Ran4PEstimator::HIterator::SplitToSecondRoutionData(double* p, std::vector< util::point2d >& fids1, std::vector<int>& idx1, std::vector< util::point2d >& fids2, 
		std::vector<int>& idx2, std::vector< std::pair< int, int > >& matches, std::vector< std::pair< int, int > >& submm, cv::Mat& HH, float dist_tol)
{
	matches.clear(); submm.clear();
	
	bool donext = true;
	
// 	CvMat* H = cvCreateMat(3, 3, CV_64FC1);
// 	cvmSet(H, 0, 0, p[0]); cvmSet(H, 0, 1, p[1]); cvmSet(H, 0, 2, p[4]);
// 	cvmSet(H, 1, 0, p[2]); cvmSet(H, 1, 1, p[3]); cvmSet(H, 1, 2, p[5]); 
// 	cvmSet(H, 2, 0, 0); cvmSet(H, 2, 1, 0); cvmSet(H, 2, 2, 1); 
	
	cv::Mat H = cv::Mat(3, 3, CV_64FC1);
    H.at<double>(0, 0)=p[0]; H.at<double>(0, 1)=p[1]; H.at<double>(0, 2)=p[4];
    H.at<double>(1, 0)=p[2]; H.at<double>(1, 1)=p[3]; H.at<double>(1, 2)=p[5];
    H.at<double>(2, 0)=0; H.at<double>(2, 1)=0; H.at<double>(2, 2)=1;
// 	std::cout<<H<<std::endl;
	
	std::vector<util::point2d> fidtrans;
	r4per->ApplyTransformToAll(fids1, fidtrans, H);
	std::vector<float> dists;
	r4per->CalculateMatch(fidtrans, fids2, matches, dists, dist_tol*.75);
	if(fids2.size()-matches.size() <= 6 || fids1.size()-matches.size() <= 6){
		r4per->CalculateMatch(fidtrans, fids2, matches, dists, dist_tol);
		donext = false;
	}
// 	if(!donext){
// 		for(int i = 0; i < fidtrans.size(); i++){
// 			std::cout<<fidtrans[i].x<<" "<<fidtrans[i].y<<std::endl;
// 		}
// 		std::cout<<std::endl;
// 		
// 		for(int i = 0; i < fids2.size(); i++){
// 			std::cout<<fids2[i].x<<" "<<fids2[i].y<<std::endl;
// 		}
// 		std::cout<<std::endl;
// 	}
// 	std::cout<<"fids1.size()"<<fids1.size()<<std::endl;
	for(int i = 0; i < matches.size(); i++){
		std::pair<int, int> mp = std::make_pair(idx1[matches[i].first], idx2[matches[i].second]);
		submm.push_back(mp);
	}
// 		std::cout<<matches.size()<<std::endl;
	std::vector<util::point2d> tmpfid1; std::vector<int> tmpidx1;
	for(int i = 0; i < fids1.size(); i++){
		bool find = false;
		for(int j = 0; j < matches.size(); j++){
			if(matches[j].first == i){
				find = true;
				break;
			}
		}
		if(!find){
			tmpfid1.push_back(fids1[i]);
			tmpidx1.push_back(idx1[i]);
		}
	}
	fids1 = tmpfid1; idx1 = tmpidx1;
	std::vector<util::point2d> tmpfid2; std::vector<int> tmpidx2;
	for(int i = 0; i < fids2.size(); i++){
		bool find = false;
		for(int j = 0; j < matches.size(); j++){
			if(matches[j].second == i){
				find = true;
				break;
			}
		}
		if(!find){
			tmpfid2.push_back(fids2[i]);
			tmpidx2.push_back(idx2[i]);
		}
	}
	fids2 = tmpfid2; idx2 = tmpidx2;
	
	r4per->ApplyTransformToAll(fids1, fidtrans, H);
	r4per->CalculateMatch(fidtrans, fids2, matches, dists, dist_tol*3);
	HH = H;
	if(matches.size() < 6){
		return false;
	}
	
	if(!donext){		
		return false;
	}
	
	r4per->EstimateTransform(fids1, fids2, matches, H);
	
	r4per->ApplyTransformToAll(fids1, fidtrans, H);
	r4per->CalculateMatch(fidtrans, fids2, matches, dists, dist_tol);

	if(matches.size() < 3){
		return false;
	}
	
// 	p[0] = H.D()[0]; p[1] = H.D()[1]; p[4] = H.D()[2]; 
// 	p[2] = H.D()[3]; p[3] = H.D()[4]; p[5] = H.D()[5]; 
    
    p[0] = H.at<double>(0, 0); p[1] = H.at<double>(0, 1); p[4] = H.at<double>(0, 2);
    p[2] = H.at<double>(1, 0); p[3] = H.at<double>(1, 1); p[5] = H.at<double>(1, 2);

	return true;
}

void Ran4PEstimator::HIterator::ReEstimateTransform(std::vector<std::vector<std::pair<int, int> > >& mmset, std::vector<cv::Mat >& hset, bool new_method)
{	
	mmset.clear();
	hset.clear();
// 	std::cout<<"*****************"<<std::endl;
	double opts[LM_OPTS_SZ], info[LM_INFO_SZ];
	
	opts[0] = LM_INIT_MU; opts[1] = 1E-15; opts[2] = 1E-15; opts[3] = 1E-20;
	opts[4] = -LM_DIFF_DELTA; // relevant only if the Jacobian is approximated using finite differences; specifies forward differencing 
	
	int m = 6;
	double p[m];
// 	p[0] = cvGetReal2D(HH, 0, 0); p[1] = cvGetReal2D(HH, 0, 1); p[4] = cvGetReal2D(HH, 0, 2);
// 	p[2] = cvGetReal2D(HH, 1, 0); p[3] = cvGetReal2D(HH, 1, 1); p[5] = cvGetReal2D(HH, 1, 2);
	
	/*p[0] = HH.D()[0]; p[1] = HH.D()[1]; p[4] = HH.D()[2]; 
	p[2] = HH.D()[3]; p[3] = HH.D()[4]; p[5] = HH.D()[5];*/ 
//     std::cout<<HH<<std::endl;
    p[0] = HH.at<double>(0, 0); p[1] = HH.at<double>(0, 1); p[4] = HH.at<double>(0, 2);
    p[2] = HH.at<double>(1, 0); p[3] = HH.at<double>(1, 1); p[5] = HH.at<double>(1, 2);
// 	std::cout<<"*****************"<<std::endl;
	std::vector<std::pair<int, int> > mm = *matches;
	
	std::vector<util::point2d> fids1 = r4per->fids1;
	std::vector<int> refidx1;
	for(int i = 0; i < fids1.size(); i++){
		refidx1.push_back(i);
	}
	std::vector<util::point2d> fids2 = r4per->fids2;
	std::vector<int> refidx2;
	for(int i = 0; i < fids2.size(); i++){
		refidx2.push_back(i);
	}
// std::cout<<"*****************"<<std::endl;
	while(true){
// 		std::cout<<1<<std::endl;
		int n = 2*mm.size();
		
		double x[n];
		for(int i = 0; i < mm.size(); i++){
			x[i*2] = fids2[mm[i].second].x;
			x[i*2+1] = fids2[mm[i].second].y;
		}
		
		std::vector<util::point2d> mfids1;
		for(int i = 0; i < mm.size(); i++){
			mfids1.push_back(fids1[mm[i].first]);
		}
		
		double* covar;
		double work[LM_DIF_WORKSZ(m, n)+m*m];
		if(!work){
			fprintf(stderr, "memory allocation request failed in main()\n");
		exit(1);
		}
		covar=work+LM_DIF_WORKSZ(m, n);
		
		dlevmar_der(HIterator::lmderfun, HIterator::lmderjac, p, x, m, n, 1000, opts, info, work, covar, &mfids1);
		std::vector<std::pair<int, int> > submm;
		cv::Mat H;
		bool donext = SplitToSecondRoutionData(p, fids1, refidx1, fids2, refidx2, mm, submm, H, r4per->pts_dist_tol);
// 		std::cout<<submm.size()<<std::endl;
		mmset.push_back(submm);
// 		std::cout<<mmset.size()<<std::endl;
		hset.push_back(H);
// 		std::cout<<refidx1.size()<<"\t"<<refidx2.size()<<std::endl;
// 		std::cout<<fids1.size()<<" "<<mm.size()<<" "<<submm.size()<<std::endl;
		if(new_method)
		{
			if(!donext){
				submm.clear();
				for(int i = 0; i < refidx1.size(); i++){
					submm.push_back(std::make_pair(refidx1[i], -1));
				}
	// 			std::cout<<refidx1.size()<<std::endl;
				mmset.push_back(submm);
				submm.clear();
				for(int i = 0; i < refidx2.size(); i++){
					submm.push_back(std::make_pair(-1, refidx2[i]));
				}
				mmset.push_back(submm);
	// 			submm.clear();
	// 			for(int i = 0; i < refidx2.size(); i++){
	// 				submm.push_back(std::make_pair(refidx2[i], -1));
	// 			}
	// 			mmset.push_back(submm);

				break;
			}
		}
		else{
			if(!donext || submm.size() == 0){
			for(int i = 0; i < refidx1.size(); i++){
				submm.push_back(std::make_pair(refidx1[i], -1));
			}
			mmset.push_back(submm);
// 			submm.clear();
// 			for(int i = 0; i < refidx2.size(); i++){
// 				submm.push_back(std::make_pair(refidx2[i], -1));
// 			}
// 			mmset.push_back(submm);

			break;
		}

		}
// 		if(!donext){
// 			submm.clear();
// 			for(int i = 0; i < refidx1.size(); i++){
// 				submm.push_back(std::make_pair(refidx1[i], -1));
// 			}
// // 			std::cout<<refidx1.size()<<std::endl;
// 			mmset.push_back(submm);
// 			submm.clear();
// 			for(int i = 0; i < refidx2.size(); i++){
// 				submm.push_back(std::make_pair(-1, refidx2[i]));
// 			}
// 			mmset.push_back(submm);
// // 			submm.clear();
// // 			for(int i = 0; i < refidx2.size(); i++){
// // 				submm.push_back(std::make_pair(refidx2[i], -1));
// // 			}
// // 			mmset.push_back(submm);
//
// 			break;
// 		}
	}
// 	std::cout<<mmset.size()<<std::endl;
}

/**xf:[a b t0; c d t1]*/
static cv::Mat xform(const util::point2d* pts1, const util::point2d* pts2, int point_size)
{
    cv::Mat H=cv::Mat::zeros(3, 3, CV_64FC1);
	int m = 2*point_size, n = 6, nrhs = 1;
	
    float As[m*n];
    float bs[m];
	
	memset(As, 0, sizeof(float)*m*n);

	float* A1 = As; float* A2 = A1+m; float* A3 = A2+m; 
	float* A4 = A3+m; float* A5 = A4+m; float* A6 = A5+m;

    for(int i = 0; i < m; i += 2){
		int pt_idx = i>>1;
		A4[i+1] = A1[i] = pts1[pt_idx].x;
		A5[i+1] = A2[i] = pts1[pt_idx].y;
		A6[i+1] = A3[i] = 1;
		bs[i] = pts2[pt_idx].x; bs[i+1] = pts2[pt_idx].y;
	}
	
	int lda = m, ldb = m;
	int jpvt[n]; memset(jpvt, 0, sizeof(int)*n);
	float rcond = -1.0;
	int rank; /* Output */
	int lwork = -1;
	float test;
	int info;

	/* Query to find a good size for the work array */
	sgelsy_(&m, &n, &nrhs, As, &lda, bs, &ldb, jpvt, &rcond, &rank, &test, &lwork, &info);

	lwork = (int)test;
	float work[lwork];

	/* Make the FORTRAN call */
	sgelsy_(&m, &n, &nrhs, As, &lda, bs, &ldb, jpvt, &rcond, &rank, work, &lwork, &info);

	/* Go from column- to row-major */
// 	for(int i = 0; i < n; i++)
// 		xf[i] = bs[i];
//     for(int i = 0; i < n; i++)
//     {
// //         std::cout<<"bs"<<i<<":"<<bs[i]<<std::endl;
//     }
    H.at<double>(0, 0)=bs[0]; H.at<double>(0, 1)=bs[1]; H.at<double>(0, 2)=bs[2];
    H.at<double>(1, 0)=bs[3]; H.at<double>(1, 1)=bs[4]; H.at<double>(1, 2)=bs[5];
    H.at<double>(2, 0)=0; H.at<double>(2, 1)=0; H.at<double>(2, 2)=1;
	if (info != 0)
		printf("Error [%d] in call to sgelsy\n", info);
//     std::cout<<H<<std::endl;
    return H;
}

static void xform_pt(const util::point2d& fid, const cv::Mat& H, util::point2d& nfid){
	const double* xf = (double*)(H.data);
	nfid.x = fid.x*xf[0]+fid.y*xf[1]+xf[2];
	nfid.y = fid.x*xf[3]+fid.y*xf[4]+xf[5];
}

void Ran4PEstimator::EstimateTransform(const std::vector< util::point2d >& fids1, 
											 const std::vector< util::point2d >& fids2, const QuadPair& qpair, cv::Mat& H)		//!WARNING
{
	util::point2d pts[3];
	util::point2d mpts[3];
	for(int i = 0; i < 3; i++){
		memcpy(pts[i].v, fids1[qpair.first[i]].v, sizeof(float)*2);
		memcpy(mpts[i].v, fids2[qpair.second[i]].v, sizeof(float)*2);
	}

	H=xform(pts, mpts, 3);
//     	std::cout<<H<<std::endl;
}

void Ran4PEstimator::EstimateTransform(const std::vector< util::point2d >& fids1, const std::vector< util::point2d >& fids2, 
									   const std::vector< std::pair< int, int > >& matches, cv::Mat& H)
{
	util::point2d pts[matches.size()];
	util::point2d mpts[matches.size()];
	for(int i = matches.size(); i--; ){
		memcpy(pts[i].v, fids1[matches[i].first].v, sizeof(float)*2);
		memcpy(mpts[i].v, fids2[matches[i].second].v, sizeof(float)*2);
	}
	
	H=xform(pts, mpts, matches.size());
}

void Ran4PEstimator::ApplyTransformToAll(const std::vector< util::point2d >& fids, std::vector< util::point2d >& fidtrans, const cv::Mat& H)
{
	fidtrans.resize(fids.size());
	for(int i = fids.size(); i--; ){
		xform_pt(fids[i], H, fidtrans[i]);
	}
}

void Ran4PEstimator::CalculateMatch(const std::vector< util::point2d >& fids1, const std::vector< util::point2d >& fids2,
									   std::vector< std::pair<int, int> >& rawmatches, std::vector<float>& dists, double dist_err_tol)
{	
	if(!fids2.size()){
		return;
	}
	
	SearchTree ref(fids2);
	int* idices = ref.IndicesArray();
	double* distances = ref.DistanceArray();
	float dist_err_thre = dist_err_tol*dist_err_tol;
	
	rawmatches.clear();
	dists.clear();
	
	for(int i = fids1.size(); i--;){
		ref.GetKNearestPoints(fids1[i], 1, idices, distances);
		if(distances[0] < dist_err_thre){
			rawmatches.push_back(std::make_pair(i, idices[0]));
			dists.push_back(distances[0]);
		}
	}
}

void Ran4PEstimator::PrepareLocalLimitedPointPairs(const std::vector< util::point2d >& fids, std::vector< Ran4PEstimator::PointPair >& ppairs)
{
	GeneratePointPairs(fids, ppairs);
	
	if(fids.size() < 42){
		return;
	}
	
	//Get the average of point-to-point minimun distance
	SearchTree distr(fids);
	int* idices = distr.IndicesArray();
	double* distances = distr.DistanceArray();
	
	float avg_dist_mean = 0;
	for(int i = fids.size(); i--;){
		distr.GetKNearestPoints(fids[i], 2, idices, distances);
		avg_dist_mean += sqrt(distances[1]);
	}
	avg_dist_mean = avg_dist_mean/fids.size();
	
	/** logic to constraint search distance **/
	float dist_mean, dist_std, min_dist, max_dist;
	PairDistanceMeanStd(ppairs, &dist_mean, &dist_std);
	min_dist = ERROR_TOL_MIN_DIST;
	
	if(avg_dist_mean < 0.5*ERROR_TOL_MIN_DIST){
		avg_dist_mean = 0.5*ERROR_TOL_MIN_DIST;
	}
	if(min_dist < avg_dist_mean){
		min_dist = avg_dist_mean;
	}
	
	max_dist = dist_mean+2*dist_std;
	
	if(max_dist > 3.5*avg_dist_mean+min_dist){			//Limit the search space so that there will not be too much combinations.
		max_dist = 3.5*avg_dist_mean+min_dist;
	}
	
	FilterPointPairs(max_dist, min_dist, ppairs);
}


void Ran4PEstimator::RansacMatch(double dist_err_tol, std::vector< std::pair< int, int > >& matches_output, 
								 std::vector<std::pair<int, int> >& ambiguity, cv::Mat& HH, bool eigenlimit)
{	
	if((int)fids1.size() < MIN_PT_NUM || (int)fids2.size() < MIN_PT_NUM){
		std::cout<<"Insufficient number of fiducial markers"<<std::endl;
		return;
	}

	PrepareLocalLimitedPointPairs(fids1, ppairs1);
	PrepareLocalLimitedPointPairs(fids2, ppairs2);
	
	pts_dist_tol = dist_err_tol;
	double squared_dist_err_tol = 0.7225*dist_err_tol*dist_err_tol;
	
	if((int)fids1.size() < MIN_PT_NUM*3 || (int)fids2.size() < MIN_PT_NUM*3){
		min_baseline = DBL_MIN;
	}
	
	int count = 0;
	float coverage = RAN_COVERAGE;
	double p = pow(1.0 - pow(coverage, 3), count);
	double error_tol = ERROR_TOL;
	
	int max_cover_num = -999;
	QuadPair max_pair;
// 	std::cout<<fids2.size()<<std::endl;

	while((p = pow(1.0 - pow(coverage, 3), count++)) > error_tol){
		fpcs::Quadrilateral quad;
		float invariant1, invariant2;
		PickRandomQuadrilateral(ppairs1, quad, &invariant1, &invariant2);

// 		std::cout<<"ref: "<<invariant1<<"\t"<<invariant2<<std::endl;

		std::vector<QuadPair> qpairs;
		SearchQuadrilateralPairs(quad, invariant1, invariant2, AREA_DIFF_TOL, pts_dist_tol*0.5, qpairs);
		for(int i = 0; i < qpairs.size(); i++){
			cv::Mat H;
			std::vector<util::point2d> fidtrans;
			std::vector<std::pair<int, int> > matches;
			std::vector<float> dists;

			EstimateTransform(fids1, fids2, qpairs[i], H);

			if(eigenlimit){											//security for H (rotation limited in 30 degree)
				float eval1, eval2;
				EigenValue2X2(H, &eval1, &eval2);
// 				if(eval1 < 0.5 || eval2 < 0.5 || eval1 > 2 || eval2 > 2){
				if(eval1 < 0.8 || eval2 < 0.8 || eval1 > 1.25 || eval2 > 1.25){
					continue;
				}
			}

			ApplyTransformToAll(fids1, fidtrans, H);
			CalculateMatch(fidtrans, fids2, matches, dists, 3*dist_err_tol);
			
			std::vector<std::pair<int, int> > varified_matches, ambiguous_matches;
			
			if(matches.size() > RAN_COVERAGE*fids1.size() && matches.size() >= 4){
				EstimateTransform(fids1, fids2, matches, H);
				ApplyTransformToAll(fids1, fidtrans, H);
				CalculateMatch(fidtrans, fids2, matches, dists, 5*dist_err_tol);
				//CalculateMatch(fidtrans, fids2, matches, dists, 3*dist_err_tol);
				
				varified_matches.reserve(matches.size());
				
				for(int i = dists.size(); i--;){
					if(dists[i] < squared_dist_err_tol){
						varified_matches.push_back(matches[i]);
					}
					else{
						ambiguous_matches.push_back(matches[i]);
					}
				}
			}
			
			if((int)varified_matches.size() > max_cover_num){
				max_cover_num = varified_matches.size();
				matches_output = varified_matches;
				ambiguity = ambiguous_matches;
				
				HH = H;
				max_pair = qpairs[i];
				coverage = ((double)max_cover_num)/fids2.size()*COVERAGE_REFRESH_RATIO;
				
// 				std::cout<<"MMATCH: "<<max_cover_num<<" "<<coverage<<std::endl;
// 				HH.Print(std::cout);
			}
		}
	}
}

void Ran4PEstimator::EigenValue2X2(const cv::Mat& H, float* eval1, float* eval2)
{
	cv::Mat A = cv::Mat::zeros(2, 2, CV_64FC1);
    A.at<double>(0, 0) = H.at<double>(0,0);
    A.at<double>(0, 1) = H.at<double>(0,1);
    A.at<double>(1, 0) = H.at<double>(1,0);
    A.at<double>(1, 1) = H.at<double>(1,1);
    
    
// 	cvmSet(A, 0, 0, cvGetReal2D(H,0,0));
// 	cvmSet(A, 0, 1, cvGetReal2D(H,0,1));
// 	cvmSet(A, 1, 0, cvGetReal2D(H,1,0));
// 	cvmSet(A, 1, 1, cvGetReal2D(H,1,1));
    
    cv::Mat evec=cv::Mat::zeros(2,2,CV_64FC1);
    cv::Mat eval=cv::Mat::zeros(2,1,CV_64FC1);
	
// 	CvMat* evec  = cvCreateMat(2,2,CV_64FC1);
// 	CvMat* eval  = cvCreateMat(2,1,CV_64FC1);
// 	cvZero(evec);
// 	cvZero(eval);
    cv::eigen(A, eval, evec);
// 	cvEigenVV(A, evec, eval, DBL_EPSILON, -1, -1);
	
    *eval1 = eval.at<double>(0,0);
    *eval2 = eval.at<double>(1,0);
// 	*eval1 = cvGetReal2D(eval,0,0);
// 	*eval2 = cvGetReal2D(eval,1,0);

}

void Ran4PEstimator::FindUnmatchedAmbiguity(const std::vector<std::vector<std::pair<int, int> > >& mmset, const std::vector<cv::Mat >& hset, double search_radius)
{
	const std::vector<std::pair<int, int> >& umidx1 = mmset[mmset.size()-2];
	const std::vector<std::pair<int, int> >& umidx2 = mmset[mmset.size()-1];
// 	std::cout<<mmset[mmset.size()-1].size()<<std::endl;
	std::vector<util::point2d> unfids1(umidx1.size()), unfids2(umidx2.size());
	
	for(int i = 0; i < umidx1.size(); i++){
		unfids1[i] = fids1[umidx1[i].first];
	}
	for(int i = 0; i < umidx2.size(); i++){
		unfids2[i] = fids2[umidx2[i].second];
       // std::cout<<"x:  "<<unfids2[i].x<<"y:  "<<unfids2[i].y<<std::endl;
	}
// std::cout<<"*****************"<<std::endl;
	SearchTree ref(unfids2);
	int* indices = ref.IndicesArray();
	double* distances = ref.DistanceArray();
	int num;
	float squared_radius = search_radius*search_radius;

	bool mask2[unfids2.size()];
	memset(mask2, 0, sizeof(bool)*unfids2.size());

	for(int j = 0; j < hset.size(); j++){
		std::vector<util::point2d> tfids;
		ApplyTransformToAll(unfids1, tfids, hset[j]);

		for(int i = tfids.size(); i--;){
			ref.GetRadiusNearestPoints(tfids[i], squared_radius, &num, indices, distances);
// 			ref.GetKNearestPoints(tfids[i], num, indices, distances);
			if(num > 0){
				unmatched_fids1.push_back(unfids1[i]);
			}
			for(int k = num; k--;){
				mask2[indices[k]] = true;
			}
// 			mask2[indices[0]] = true;
		}
	}
	for(int i = 0; i < unfids2.size(); i++){
		if(mask2[i]){
			unmatched_fids2.push_back(unfids2[i]);
		}
	}
// 		std::cout<<"??????????"<<std::endl;
}

void Ran4PEstimator::AffineTransformEstimation(double dist_err_tol, std::vector<std::pair<util::point2d, util::point2d> >& matchvec, std::vector<cv::Mat >& hset, bool test, bool eigenlimit)
{
	std::vector<std::pair<int, int> > matches, ambiguity;

	cv::Mat HH;
	RansacMatch(dist_err_tol, matches, ambiguity, HH, eigenlimit);
// 	std::cout<<"HH"<<std::endl;
	HIterator::SetInitialInput(this, HH, &matches);
//     std::cout<<HH<<std::endl;
	std::vector<std::vector<std::pair<int, int> > > mmset;
	HIterator::ReEstimateTransform(mmset, hset, true);
// std::cout<<"*****************"<<std::endl;
	matches.clear();
// 	std::cout<<mmset[mmset.size()-2].size()<<std::endl;
	for(int i = 0; i < mmset.size()-2; i++){		// the last two subsets are the unmatched fiducials in fids1 and fids2
		std::vector<std::pair<int, int> >& mm = mmset[i];
		for(int j = 0; j < mm.size(); j++){
// 			std::cout<<mm[j].first<<std::endl;
			matchvec.push_back(std::make_pair(fids1[mm[j].first], fids2[mm[j].second]));
             //std::cout<<"first:"<<mm[j].first<<"second:"<<mm[j].second<<std::endl;
		}
	}

// 	for(int i=0;i,mmset.size()-2;i++)
//     {
//         std::vector<std::pair<int, int> >& nn = mmset[i];
//         for(int j;j<nn.size();j++)
//         {
//             std::cout<<"first:"<<nn[j].first<<"second:"<<nn[j].second<<std::endl;
//         }
//     }
	FindUnmatchedAmbiguity(mmset, hset, 3.5*pts_dist_tol);
// 	std::cout<<"***************"<<std::endl;
	if(unmatched_fids1.size() < 1 || unmatched_fids2.size() < 1){	//the match is perfect
		return;
	}

#define ANK_NUM			4	

	for(int i = 1; i <= ANK_NUM; i++){		//for the stability of drift estimation 
		unmatched_fids1.push_back(matchvec[i].first);
		unmatched_fids2.push_back(matchvec[i].second);
	}
}

std::vector< util::point2d >& Ran4PEstimator::UmFids1()
{
	return unmatched_fids1;
}

std::vector< util::point2d >& Ran4PEstimator::UmFids2()
{
	return unmatched_fids2;
}

void Ran4PEstimator::DistributionCorrection(double dist_err_tol, std::vector< std::pair< util::point2d, util::point2d > >& matchvec)
{
	SearchTree ref1(fids1);
	SearchTree ref2(fids2);
	int* idices = ref1.IndicesArray();
	double* distances = ref1.DistanceArray();
	float dist_err_thre = dist_err_tol*dist_err_tol;
	
	int count = 0;
	
	for(std::vector< std::pair< util::point2d, util::point2d > >::iterator itr = matchvec.begin(); itr != matchvec.end();){
		ref1.GetKNearestPoints((*itr).first, 4, idices, distances);
		if(distances[3] < dist_err_thre){
			itr = matchvec.erase(itr);
			count++;
		}
		else{
			itr++;
		}
	}
	
	for(std::vector< std::pair< util::point2d, util::point2d > >::iterator itr = matchvec.begin(); itr != matchvec.end();){
		ref2.GetKNearestPoints((*itr).second, 4, idices, distances);
		if(distances[3] < dist_err_thre){
			itr = matchvec.erase(itr);
			count++;
		}
		else{
			itr++;
		}
	}
	std::cout<<count<<" removed."<<std::endl;
}

void Ran4PEstimator::GenerateQuadrilateralVec(std::vector<fpcs::Quadrilateral>& quadvec, const std::vector<util::point2d>& fids,
											  std::vector<std::pair<float, float> >& invavec, std::vector<std::pair<float, float> >& baselinevec)
{
	for(int i = 0; i < fids.size(); i++){
		for(int j = i+1; j < fids.size(); j++){
			for(int k = j+1; k < fids.size(); k++){
				for(int l = k+1; l < fids.size(); l++){
					fpcs::Quadrilateral quad(i,j,k,l);
					float invariant1, invariant2;
					float baseline1, baseline2;
					int index[4];
					if(fpcs::QuadrilateralInvariant(fids[quad[0]], fids[quad[1]], fids[quad[2]], fids[quad[3]], index, &invariant1, &invariant2, &baseline1, &baseline2)
					&&invariant1>=quad_svd_safe_tol&&invariant1<=1-quad_svd_safe_tol&&invariant2>=quad_svd_safe_tol&&invariant2<=1-quad_svd_safe_tol){
						index[0] = quad[index[0]];
						index[1] = quad[index[1]];
						index[2] = quad[index[2]];
						index[3] = quad[index[3]];
						quad[0] = index[0];
						quad[1] = index[1];
						quad[2] = index[2];
						quad[3] = index[3];
						quadvec.push_back(quad);
						invavec.push_back(std::make_pair(invariant1, invariant2));
						baselinevec.push_back(std::make_pair(baseline1, baseline2));
					}
				}
			}
		}
	}
}

void Ran4PEstimator::RefineQuadrilateralVec(std::vector< fpcs::Quadrilateral >& quadvec, std::vector< std::pair< float, float > >& invavec,
											std::vector< std::pair< float, float > >& baselinevec, double* graph_diameter, double* min_baseline)
{
	*graph_diameter = 0;
	for(int i = 0; i < baselinevec.size(); i++){
		if(*graph_diameter < baselinevec[i].first){
			*graph_diameter = baselinevec[i].first;
		}
		if(*graph_diameter < baselinevec[i].second){
			*graph_diameter = baselinevec[i].second;
		}
	}
	*min_baseline = *graph_diameter*graph_baseline_safe_tol;
	std::vector<fpcs::Quadrilateral> tquadvec;
	std::vector< std::pair< float, float > > tinvavec;
	std::vector< std::pair< float, float > > tbaselinevec;

	for(int i = 0; i < quadvec.size(); i++){
		if(baselinevec[i].first >= *min_baseline && baselinevec[i].second >= *min_baseline){
			tquadvec.push_back(quadvec[i]);
			tinvavec.push_back(invavec[i]);
			tbaselinevec.push_back(baselinevec[i]);
		}
	}
	quadvec = tquadvec;
	invavec = tinvavec;
	baselinevec = tbaselinevec;
}

void Ran4PEstimator::SetDefaultSearchParams(int __nk, double __eps, double __quad_svd_safe_tol, double __baseline_safe_tol)
{
	nk = __nk;
	eps = __eps;
	quad_svd_safe_tol = __quad_svd_safe_tol;
	graph_baseline_safe_tol = __baseline_safe_tol;
	if(nnIdx){
		delete [] nnIdx;
		delete [] dists;
		nnIdx = new ANNidx[nk];						// allocate near neigh indices
		dists = new ANNdist[nk];					// allocate near neighbor dists
	}
}

void Ran4PEstimator::CreateReferenceTree(const std::vector< std::pair< float, float > >& invavec)
{
	nPts = invavec.size();
	queryPt = annAllocPt(dim);
	dataPts = annAllocPts(nPts, dim);			// allocate data points
	if(nnIdx == NULL && dists == NULL){
		nnIdx = new ANNidx[nk];						// allocate near neigh indices
		dists = new ANNdist[nk];						// allocate near neighbor dists
	}
    for(int i = 0; i < invavec.size(); i++){
		dataPts[i][0] = invavec[i].first;
		dataPts[i][1] = invavec[i].second;
	}
	kdTree = new ANNkd_tree(dataPts, nPts, dim);
	if(nk > nPts){
		nk = nPts-1;
	}
}

void Ran4PEstimator::GenerateRandomQuadrilateral(fpcs::Quadrilateral& quad, const std::vector< util::point2d >& fids,
												 float* invariant1, float* invariant2, float* baseline1, float* baseline2)
{
	bool sucess = false;
	while(!sucess){
		fpcs::GenerateRandomQuadIndex(fids.size(), &quad);

		int index[4];
		if(fpcs::QuadrilateralInvariant(fids[quad[0]], fids[quad[1]], fids[quad[2]], fids[quad[3]], index, invariant1, invariant2, baseline1, baseline2)
		&& *invariant1 >= quad_svd_safe_tol && *invariant1 <= 1-quad_svd_safe_tol && *invariant2 >= quad_svd_safe_tol && *invariant2 <= 1-quad_svd_safe_tol
		&& *baseline1 >= min_baseline && *baseline2 >= min_baseline){
			index[0] = quad[index[0]];
			index[1] = quad[index[1]];
			index[2] = quad[index[2]];
			index[3] = quad[index[3]];
			quad[0] = index[0];
			quad[1] = index[1];
			quad[2] = index[2];
			quad[3] = index[3];
			sucess = true;
		}
	}
}

void Ran4PEstimator::SearchQuadrilateralPairs_Global(const fpcs::Quadrilateral& quad, float invariant1, float invariant2, std::vector<QuadPair>& pairs)
{
	pairs.clear();

	fpcs::Quadrilateral quads[4];
	float invariant1s[4];
	float invariant2s[4];

	quads[0][0] = quad[0];quads[0][1] = quad[1];quads[0][2] = quad[2];quads[0][3] = quad[3];
	invariant1s[0] = invariant1; invariant2s[0] = invariant2;

	quads[1][0] = quad[1];quads[1][1] = quad[0];quads[1][2] = quad[2];quads[1][3] = quad[3];
	invariant1s[1] = 1-invariant1; invariant2s[1] = invariant2;

	quads[2][0] = quad[0];quads[2][1] = quad[1];quads[2][2] = quad[3];quads[2][3] = quad[2];
	invariant1s[2] = invariant1; invariant2s[2] = 1-invariant2;

	quads[3][0] = quad[1];quads[3][1] = quad[0];quads[3][2] = quad[3];quads[3][3] = quad[2];
	invariant1s[3] = 1-invariant1; invariant2s[3] = 1-invariant2;

	for(int i = 0; i < 4; i++){
		queryPt[0] = invariant1s[i];
		queryPt[1] = invariant2s[i];
		kdTree->annkSearch(queryPt, nk, nnIdx, dists, eps);
		for(int j = 0; j < nk; j++){
			QuadPair qpair;
			qpair.first = quads[i];
			qpair.second = quadvec[nnIdx[j]];
			pairs.push_back(qpair);
		}
	}
}


void Ran4PEstimator::RansacMatch_Global(double dist_err_tol, std::vector<std::pair<int, int>>& matches_output, cv::Mat& HH, bool eigenlimit)
{
	if((int)fids1.size() < MIN_PT_NUM || (int)fids2.size() < MIN_PT_NUM){
		std::cout<<"Insufficient number of fiducial markers"<<std::endl;
		return;
	}
	else if((int)fids1.size() < MIN_PT_NUM*3 || (int)fids2.size() < MIN_PT_NUM*3){
		quad_svd_safe_tol = 0;
		graph_baseline_safe_tol = 0;
	}

	pts_dist_tol = dist_err_tol;

	GenerateQuadrilateralVec(quadvec, fids2, invavec, baselinevec);
	RefineQuadrilateralVec(quadvec, invavec, baselinevec, &graph_diameter, &min_baseline);

	if((int)fids1.size() < MIN_PT_NUM*3 || (int)fids2.size() < MIN_PT_NUM*3){
		min_baseline = DBL_MIN;
	}

	SetDefaultSearchParams(KDTREE_NK, 0, quad_svd_safe_tol, graph_baseline_safe_tol);
	CreateReferenceTree(invavec);

	int count = 0;
	float coverage = RAN_COVERAGE_old;
	double p = pow(1.0 - pow(coverage, 4), count);
	double error_tol = ERROR_TOL_old;

	int max_cover_num = -999;
	HH = NULL;
	QuadPair max_pair;
	while(p = pow(1.0 - pow(coverage, 4), count++) > error_tol){
// 		std::cout<<p<<"("<<error_tol<<")"<<":"<<count<<std::endl;
		fpcs::Quadrilateral quad;
		float invariant1, invariant2;
		float baseline1, baseline2;

		GenerateRandomQuadrilateral(quad, fids1, &invariant1, &invariant2, &baseline1, &baseline2);

		std::vector<QuadPair> qpairs;
		SearchQuadrilateralPairs_Global(quad, invariant1, invariant2, qpairs);
		for(int i = 0; i < qpairs.size(); i++){
			cv::Mat H;
			std::vector<util::point2d> fidtrans;
			std::vector<std::pair<int, int> > matches;
			std::vector<float> dists;

			EstimateTransform(fids1, fids2, qpairs[i], H);

// 			if((cvGetReal2D(H,0,0) < 0.5 && cvGetReal2D(H,1,1) < 0.5) ||
// 				(cvGetReal2D(H,0,0) > 1.5 && cvGetReal2D(H,1,1) > 1.5)){		//security for H WARNING
// 				cvReleaseMat(&H);
// 				continue;
// 			}

			if(eigenlimit){											//security for H (rotation limited in 30 degree)
				float eval1, eval2;
				EigenValue2X2(H, &eval1, &eval2);
				if(eval1 < 0.8 || eval2 < 0.8 || eval1 > 1.25 || eval2 > 1.25){
					continue;
				}
			}

			ApplyTransformToAll(fids1, fidtrans, H);
			CalculateMatch(fidtrans, fids2, matches, dists, dist_err_tol*3);

// 			std::cout<<matches.size()<<std::endl;

			if(matches.size() > RAN_COVERAGE_old*fids1.size() && matches.size() >= 4){
				EstimateTransform(fids1, fids2, matches, H);
				ApplyTransformToAll(fids1, fidtrans, H);
				CalculateMatch(fidtrans, fids2, matches, dists, dist_err_tol);
			}

			if((int)matches.size() > max_cover_num){
				max_cover_num = matches.size();
				matches_output = matches;

				HH = H;
				max_pair = qpairs[i];
				coverage = ((double)max_cover_num)/fids2.size()*COVERAGE_REFRESH_RATIO_old;
			}

		}
	}

}



void Ran4PEstimator::GlobalMatch(double dist_err_tol, std::vector<std::pair<util::point2d, util::point2d>>& matchvec, std::vector<cv::Mat>& hset, bool test, bool eigenlimit)
{
	SetPointSetDistErrorTOL(dist_err_tol);
	std::vector<std::pair<int, int> > matches;
	cv::Mat HH;
	RansacMatch_Global(dist_err_tol, matches, HH, eigenlimit);

// 	std::cout<<"HH"<<std::endl;
// std::cout<<matches[0].first<<std::endl;
	if(!HH.empty())
	{
		HIterator::SetInitialInput(this, HH, &matches);
		std::vector<std::vector<std::pair<int, int> > > mmset;
		HIterator::ReEstimateTransform(mmset, hset, false);
// 		std::cout<<mmset.size()<<std::endl;
		matches.clear();

// 		for(int i = 0; i < mmset.size()-1; i++){		// the last subset is the unmatched fiducials in fids1
// 		std::vector<std::pair<int, int> >& mm = mmset[i];
// 		for(int j = 0; j < mm.size(); j++){
// 			matchvec.push_back(std::make_pair(fids1[mm[j].first], fids2[mm[j].second]));
// 		}
// 	}

		for(int i = 0; i < mmset.size()-1; i++){		// the last subset is the unmatched fiducials in fids1
			std::vector<std::pair<int, int> >& mm = mmset[i];
// 			std::cout<<mm.size()<<std::endl;
			for(int j = 0; j < mm.size(); j++){
// 				std::cout<<mm.size()<<std::endl;
				matchvec.push_back(std::make_pair(fids1[mm[j].first], fids2[mm[j].second]));
// 				std::cout<<j<<"  "<<fids2[mm[j].second].x<<std::endl;
			}
		}
	}
// 	std::cout<<matchvec.size()<<std::endl;

// 	HIterator::SetInitialInput(this, HH, &matches);
// 	std::vector<std::vector<std::pair<int, int> > > mmset;
// 	HIterator::ReEstimateTransform(mmset, hset);
//
// 	matches.clear();
//
// 	for(int i = 0; i < mmset.size()-1; i++){		// the last subset is the unmatched fiducials in fids1
// 		std::vector<std::pair<int, int> >& mm = mmset[i];
// 		for(int j = 0; j < mm.size(); j++){
// 			matchvec.push_back(std::make_pair(fids1[mm[j].first], fids2[mm[j].second]));
// 		}
// 	}

	if(test){

	}

}



// Ran4PEstimator::SearchTree::SearchTree(const std::vector< util::point2d >& fixed, bool alloc_buffer):nk(fixed.size()), psdim(2), pseps(0.01), psNPts(fixed.size())
// {
//     psQryPt = annAllocPt(psdim); // query point
//     psDataPts = annAllocPts(psNPts, psdim); // data points
// 	if(alloc_buffer){
// 		psNNIdx = new ANNidx[nk]; // near neighbor indices
// 		psDists = new ANNdist[nk]; // near neighbor distances
// 	}
// 	else{
// 		psNNIdx = NULL;
// 		psDists = NULL;
// 	}
// 
//     for(int i = fixed.size(); i--;){
// 		psDataPts[i][0] = fixed[i].x;
// 		psDataPts[i][1] = fixed[i].y;
// 	}
// 	
// 	psKdTree = new ANNkd_tree(psDataPts, psNPts, psdim);
// }
// 
// double* Ran4PEstimator::SearchTree::DistanceArray()
// {
// 	return psDists;
// }
// 
// int* Ran4PEstimator::SearchTree::IndicesArray()
// {
// 	return psNNIdx;
// }
// 
// void Ran4PEstimator::SearchTree::GetKNearestPoints(const util::point2d& quary, int k, int* indices, double* squared_distances)
// {
// 	psQryPt[0] = quary.x;
// 	psQryPt[1] = quary.y;
// 	psKdTree->annkSearch(psQryPt, k, indices, squared_distances, pseps);
// }
// 
// void Ran4PEstimator::SearchTree::GetRadiusNearestPoints(const util::point2d& quary, double radius, int* k, int* indices, double* squared_distances)
// {
// 	psQryPt[0] = quary.x;
// 	psQryPt[1] = quary.y;
// 	*k = psKdTree->annkFRSearch(psQryPt, radius, psNPts, indices, squared_distances, pseps);
// }
// 
// Ran4PEstimator::SearchTree::~SearchTree()
// {
// 	if(psNNIdx) delete [] psNNIdx;
//     if(psDists) delete [] psDists;
//     if(psKdTree) delete psKdTree;
//     if(psDataPts) annDeallocPts(psDataPts);
//     if(psQryPt) annDeallocPt(psQryPt);
// }


