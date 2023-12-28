#include "cpdestimator.h"
#include "util/exception.h"
#include "matrix/matrix.h"
#include "cmath"

static void PrintPointSet(const std::vector<util::point2d>& pts)
{
	for(int i = 0; i < pts.size(); i++){
		std::cout<<pts[i].x<<" "<<pts[i].y<<std::endl;
	}
	std::cout<<std::endl;
}

void CPDEstimator::Normalizer::CalculateCenter(const std::vector< util::point2d >& pts, util::point2d& center)
{
	center.x = center.y = 0;		//calcuate the center of point set
	for(size_t i = pts.size(); i--;){  //for(size_t i = 0; i < pts.size(); i++){
		center.x += pts[i].x; center.y += pts[i].y;
	}
	
	center.x /= pts.size();
	center.y /= pts.size();
}

void CPDEstimator::Normalizer::CalculateScale(const std::vector< util::point2d >& pts, float* scale)
{
	//scale
	*scale = 0;
	for(size_t i = pts.size(); i--;){  //for(size_t i = 0; i < pts.size(); i++){
		*scale += pts[i].x*pts[i].x+pts[i].y*pts[i].y;
	}
	
	*scale /= pts.size();
	*scale = sqrt(*scale);
}

void CPDEstimator::Normalizer::MovePointSetToCenter(const util::point2d& center, std::vector< util::point2d >& pts)
{
	//move the point set to center
	for(size_t i = pts.size(); i--;){  //for(size_t i = 0; i < pts.size(); i++){		
		pts[i].x += center.x; pts[i].y += center.y;
	}
}

void CPDEstimator::Normalizer::ScalePointSet(const float scale, std::vector< util::point2d >& pts)
{
	for(size_t i = pts.size(); i--;){  //for(size_t i = 0; i < pts.size(); i++){		
		pts[i].x /= scale; pts[i].y /= scale;
	}
}

void CPDEstimator::Normalizer::GlobalNormalize(std::vector< util::point2d >& fixed, std::vector< util::point2d >& moving, 
									util::point2d& fixed_reverse_center, util::point2d& moving_reverse_center, float* scale)
{
// 	CalculateCenter(fixed, fixed_reverse_center);
// 	CalculateCenter(moving, moving_reverse_center);
// 	
// 	MovePointSetToCenter(util::point2d(-fixed_reverse_center.x, -fixed_reverse_center.y), fixed);
// 	MovePointSetToCenter(util::point2d(-moving_reverse_center.x, -moving_reverse_center.y), moving);
	
	memset(fixed_reverse_center.v, 0, sizeof(float)*2);
	memset(moving_reverse_center.v, 0, sizeof(float)*2);
	
	float fixed_scale, moving_scale;
	CalculateScale(fixed, &fixed_scale);
	CalculateScale(moving, &moving_scale);
	*scale = max(fixed_scale, moving_scale);
	
	ScalePointSet(*scale, fixed);
	ScalePointSet(*scale, moving);
}

void CPDEstimator::Normalizer::ReverseNormalization(const util::point2d& reverse_center, float scale, std::vector< util::point2d >& points)
{
	ScalePointSet(1/scale, points);
	MovePointSetToCenter(reverse_center, points);
}


double CPDEstimator::DefaultSigma2(const std::vector< util::point2d >& fixed, const std::vector< util::point2d >& moving)
{
	double tf = 0, tm = 0, sumfx = 0, sumfy = 0, summx = 0, summy = 0;
	int dim = 2;
	
	for(size_t i = fixed.size(); i--;){		//for(size_t i = 0; i < fixed.size(); i++){
		tf += fixed[i].x*fixed[i].x+fixed[i].y*fixed[i].y;
		sumfx += fixed[i].x; sumfy += fixed[i].y;
	}
	
	for(size_t i = moving.size(); i--;){		//for(size_t i = 0; i < moving.size(); i++){
		tm += moving[i].x*moving[i].x+moving[i].y*moving[i].y;
		summx += moving[i].x; summy += moving[i].y; 
	}
	
	return (moving.size()*tf+fixed.size()*tm-2*(sumfx*summx+sumfy*summy))/(fixed.size()*moving.size()*dim);
}

inline static float PT2DISTANCE(const util::point2d& a, const util::point2d& b)
{
	float deltax = ((a).x-(b).x);
	float deltay = ((a).y-(b).y);
	return deltax*deltax+deltay*deltay;
}

/** Given fixed PointSet (X) = the data points; moving PointSet (Y) = the GMM centroids;
 * 	X has N elements and Y has M elements;
 *  
 *  The GMM probability density function is:
 * 
 *  p(x)=sum_{m=1}^{M+1}{P(m)p(X|m)}
 *  p(x|m)=(2*pi*sigma^2 )^{-D/2}*exp{-(x-y_{m})^2/(2*sigma^2)}, D is the dimension
 *  P(m) respect the choice to choose the $m$th centroid; P(m)=M^{-1}
 * 
 * 	additional distribution p(x|M+1)=N^{-1} is added to represent noise and outliers
 *  priori(先验) for dataset and noise is defined as uniform distribution
 *  given weight w, we have:
 * 
 *  p(x)=w/N+(1-w)*sum_{m}^{M}{1/M*p(x|m)}
 * 
 *  use theta, sigma^2 denote the model, the MLE aim is:
 *  
 *  E(theta,sigma^2)=pow_{n=1}^{N}{p(x)}
 * 
 *  the log version of MLE aim is:
 *  
 *  E(theta,sigma^2)=sum_{n=1}^{N}{log{sum_{m=1}^{M+1}{P(m)p(X|m)}}}
 * 
 *  the correspondence probability between y_m and x_n:
 * 
 *  P(m|x_n)=P(m)p(x_n|m)/p(x_n)
 * 
 *  as EM method, our Estimation period: P^{old}(m|x_n)=P(m)p(x_n|m)/p(x_n)
 *  
 *  Q=-sum_{n=1}^{N}{sum_{m=1}^{M+1}{P^{old}(m|x_n)log(P^{new}(m)p^{new}(x_n|m))}}
 * 
 *  Q(theta,sigma^2)=(2*sigma^2)^{-1}*sum_{n=1}^{N}{sum_{m=1}^{M}{P^{old}(m|x_n)*(x_n-T(y_m,theta))^2}}+(N_p*D)/2*log(sigma^2)
 * 
 *  N_p=sum_{n=1}^{N}{sum_{m=1}^{M}{P^{old}(m|x_n)}} <= N
 * 
 *  P^{old}(m|x_n)=exp(-2^{-1}*((x_n-T(y_m,theta^{old}))/sigma)^2)/(sum_{k=1}^M{exp(-2^{-1}*((x_n-T(y_m,theta^{old}))/sigma)^2)}+c)
 * 
 *  c=(2*pi*sigma^2)^{D/2}*{w/(1-w)}*M/N
 * 
 * */
void CPDEstimator::ComputeDirectGaussTransform(const std::vector< util::point2d >& fixed, 
										  const std::vector< util::point2d >& moving, double sigma2, double outliers, Probabilities& probabilities)
{
	double ksig = -2.0 * sigma2;		//sigma2 <--> sigma^2
	double ksig_1 = 1/ksig;
    int dim = 2;				//dimension
    double c = (outliers*moving.size()*std::pow(-ksig*M_PI, 0.5*dim))/((1-outliers)*fixed.size());		// c = (w*M*(2*pi*sigma^2)^{D/2})/((1-w)*N)
    
    std::vector<double> p(moving.size(), 0);
	
    std::vector<double>& p1 =  *(new(&(probabilities.p1))std::vector<double>(moving.size(), 0));
    std::vector<double>& pt1 = *(new(&(probabilities.pt1))std::vector<double>(fixed.size(), 0));
	std::vector<util::point2d>& px = *(new(&(probabilities.px))std::vector<util::point2d>(moving.size(), util::point2d(0.0f, 0.0f)));
    double& l = probabilities.l = 0.0;
	
    for(size_t i = fixed.size(); i--;){			//for(size_t i = 0; i < fixed.size(); i++){		// n=1...N
        double sp = 0;

        for(size_t j = moving.size(); j--;){			//for(size_t j = 0; j < moving.size(); j++){		// m=1...M
            double razn = PT2DISTANCE(fixed[i], moving[j]); 	// (x_n-y_m)^2
            p[j] = std::exp(razn*ksig_1);			// p_{n=i,m=j}=exp{-(x_{n}-y_{m})^2/(2*sigma^2)}
            sp += p[j];
        }
        sp += c;		// sp_{n=i}=sum_{m=1}^{M}{exp{-(x_{n}-y_{m})^2/(2*sigma^2)}}+c
        double sp_1 = 1/sp;
        pt1[i] = 1 - c*sp_1;		// pt1_{n=i}=(sum_{m=1}^{M}{exp{-(x_{n}-y_{m})^2/(2*sigma^2)}})/(sum_{m=1}^{M}{exp{-(x_{n}-y_{m})^2/(2*sigma^2)}}+c)
        
        for(size_t j = moving.size(); j--;){			//for(size_t j = 0; j < moving.size(); ++j){		// m=1...M
			double pjsp_1 = p[j]*sp_1;
            p1[j] += pjsp_1;			//finally, p1_{m=j}=sum_{n=1}^{N}{exp{-(x_{n}-y_{m})^2/(2*sigma^2)}/(sum_{k=1}^{M}{exp{-(x_{n}-y_{k})^2/(2*sigma^2)}}+c)}
            px[j].x += fixed[i].x*pjsp_1; px[j].y += fixed[i].y*pjsp_1;
        }
        l += -std::log(sp);
    }
    l += dim * fixed.size() * std::log(sigma2) / 2;		//l <--> Q
//     p1.clear();
//     pt1.clear();
//     px.clear();
}

void CPDEstimator::CalculateCorrespondence(const std::vector< util::point2d >& fixed, 
										  const std::vector< util::point2d >& moving, double sigma2, double outliers, std::vector<int>& correspondence)
{
	double ksig = -2.0 * sigma2;		//sigma2 <--> sigma^2
	double ksig_1 = 1/ksig;
    int dim = 2;				//dimension
    double c = (outliers*moving.size()*std::pow(-ksig*M_PI, 0.5*dim))/((1-outliers)*fixed.size());		// c = (w*M*(2*pi*sigma^2)^{D/2})/((1-w)*N)
    
    std::vector<double> p(moving.size(), 0);
	std::vector<double> p1_max(moving.size(), 0);
	
    new(&correspondence)std::vector<int>(moving.size());
	
    for(size_t i = fixed.size(); i--;){			//for(size_t i = 0; i < fixed.size(); i++){		// n=1...N
        double sp = 0;

        for(size_t j = moving.size(); j--;){			//for(size_t j = 0; j < moving.size(); j++){		// m=1...M
            double razn = PT2DISTANCE(fixed[i], moving[j]); 	// (x_n-y_m)^2
            p[j] = std::exp(razn*ksig_1);			// p_{n=i,m=j}=exp{-(x_{n}-y_{m})^2/(2*sigma^2)}
            sp += p[j];
        }
        sp += c;		// sp_{n=i}=sum_{m=1}^{M}{exp{-(x_{n}-y_{m})^2/(2*sigma^2)}}+c
        double sp_1 = 1/sp;
        
        for(size_t j = moving.size(); j--;){			//for(size_t j = 0; j < moving.size(); ++j){		// m=1...M
			double pjsp_1 = p[j]*sp_1;
			
            if(pjsp_1 > p1_max[j]){	//p[j]/sp<-->P(m=j|x_{n=i}); p1_max[j]<-->Maximum possibility 
                correspondence[j] = i;
                p1_max[j] = pjsp_1;
            }
        }
    }
}

double CPDEstimator::ComputeTreeGaussCorrelation(const std::vector< util::point2d >& fixed, const std::vector< util::point2d >& moving, double sigma2)
{
	double ksig = -2.0 * sigma2;		//sigma2 <--> sigma^2
	double thre = log(1.0 / ERROR_TOL)*sigma2;			//3*sigma threshold
	double ksig_1 = 1/ksig;
	double sp = 0;
	int k;
	int* indices = ftree.IndicesArray();
	double* dists = ftree.DistanceArray();
	for(size_t j = moving.size(); j--;){
		ftree.GetRadiusNearestPoints(moving[j], thre, &k, indices, dists);
		for(int i = k; i--;){
			sp += std::exp(dists[i]*ksig_1);
		}
	}
	return sp;
}

double CPDEstimator::ComputeDirectGaussCorrelation(const std::vector< util::point2d >& fixed, const std::vector< util::point2d >& moving, double sigma2)
{
	double ksig = -2.0 * sigma2;		//sigma2 <--> sigma^2
	double thre = 9 * sigma2;			//3*sigma threshold
	double ksig_1 = 1/ksig;
	double sp = 0;
	
	for(size_t i = fixed.size(); i--;){
		for(size_t j = moving.size(); j--;){
            double razn = PT2DISTANCE(fixed[i], moving[j]); 	// (x_n-y_m)^2
			if(razn > thre){
				continue;
			}
			sp += std::exp(razn*ksig_1);
		}
	}
	return sp;
}

void CPDEstimator::ComputeFastGaussTransform(const std::vector< util::point2d >& fixed, 
											 const std::vector< util::point2d >& moving, double sigma2, double outliers, Probabilities& prob)
{
	//sigma2 <--> sigma^2
	double ksig = -2.0 * sigma2;
    int dim = 2;		//dimension
    
    std::vector<double> a;
    
    ifgt::IFGaussianTransform ifgtm2f(moving, sqrt(2*sigma2), ERROR_TOL);
	ifgtm2f.CalculateTransformWith(fixed, a);
    double c = (outliers*moving.size()*std::pow(-ksig*M_PI, 0.5*dim))/((1-outliers)*fixed.size());		// c = (w*M*(2*pi*sigma^2)^{D/2})/((1-w)*N)
	
	new (&(prob.pt1))std::vector<double>(fixed.size());
	prob.l = 0;
	for(int i = a.size(); i--;){				//for(int i = 0; i < a.size(); i++){
		prob.pt1[i] = a[i]/(a[i]+c);
		prob.l += -std::log(a[i]+c);
		a[i] = 1/(a[i]+c);
	}
	prob.l += dim*fixed.size()*std::log(sigma2)/2;		//l <--> Q
	
	ifgt::IFGaussianTransform ifgtf2m(fixed, sqrt(2*sigma2), ERROR_TOL);
	ifgtf2m.CalculateTransformWith(moving, a, prob.p1);
	
	{
		std::vector<double> axy(fixed.size());
		for(int i = axy.size(); i--;){				//for(int i = 0; i < axy.size(); i++){
			axy[i] = a[i]*fixed[i].x;
		}
		
		std::vector<double> tmprslt;
		new (&(prob.px))std::vector<util::point2d>(moving.size());
		ifgtf2m.CalculateTransformWith(moving, axy, tmprslt);
		for(int i = prob.px.size(); i--;){			//for(int i = 0; i < prob.px.size(); i++){
			prob.px[i].x = tmprslt[i];
		}
		
		for(int i = axy.size(); i--;){		//for(int i = 0; i < axy.size(); i++){
			axy[i] = a[i]*fixed[i].y;
		}
		
		ifgtf2m.CalculateTransformWith(moving, axy, tmprslt);
		for(int i = prob.px.size(); i--;){		//for(int i = 0; i < prob.px.size(); i++){
			prob.px[i].y = tmprslt[i];
		}
	}
}

void CPDEstimator::MaximumStepAsAffineTransform(const std::vector<util::point2d>& fixed, const std::vector<util::point2d>& moving,
											   const Probabilities& probabilities, double sigma2, Params& result){
    int dim = 2;
	double np = 0;
	
	for(size_t i = probabilities.pt1.size(); i--;){	
		np += probabilities.pt1[i];
	}
	
	util::point2d mu_x(0,0), mu_y(0,0);
	
	for(size_t i = fixed.size(); i--;){	
		mu_x.x += fixed[i].x*probabilities.pt1[i]; mu_x.y += fixed[i].y*probabilities.pt1[i];
	}
	mu_x.x /= np; mu_x.y /= np;
	
	for(size_t i = moving.size(); i--;){
		mu_y.x += moving[i].x*probabilities.p1[i]; mu_y.y += moving[i].y*probabilities.p1[i];
	}
	mu_y.x /= np; mu_y.y /= np;

	double b1[4], b2[4], tmp[4];
	memset(b1, 0, sizeof(double)*4);
	memset(b2, 0, sizeof(double)*4);
	
	for(size_t i = moving.size(); i--;){
		b1[0] += probabilities.px[i].x*moving[i].x;
		b1[1] += probabilities.px[i].x*moving[i].y;
		b1[2] += probabilities.px[i].y*moving[i].x;
		b1[3] += probabilities.px[i].y*moving[i].y;
		
		b2[0] += probabilities.p1[i]*moving[i].x*moving[i].x;
		b2[1] += probabilities.p1[i]*moving[i].y*moving[i].x;
		b2[3] += probabilities.p1[i]*moving[i].y*moving[i].y;
	}
	
	b1[0] -= np*mu_x.x*mu_y.x;
	b1[1] -= np*mu_x.x*mu_y.y;
	b1[2] -= np*mu_x.y*mu_y.x;
	b1[3] -= np*mu_x.y*mu_y.y;
	
	b2[0] -= np*mu_y.x*mu_y.x;
	b2[1] -= np*mu_y.x*mu_y.y;
	b2[2] = b2[1];
	b2[3] -= np*mu_y.y*mu_y.y;
	
	double* A = &(result.p[0]);
	double* t = &(result.p[4]);
	
	matrix_invert_inplace(2, b2);
	matrix_product(2, 2, 2, 2, b1, b2, A);
	t[0] = mu_x.x - (A[0]*mu_y.x+A[1]*mu_y.y);
	t[1] = mu_x.y - (A[2]*mu_y.x+A[3]*mu_y.y);
	
	result.scale = 1.0;
	
	result.sigma2 = 0;
	for(size_t i = fixed.size(); i--;){			//for(size_t i = 0; i < fixed.size(); i++){
		result.sigma2 += (fixed[i].x*fixed[i].x+fixed[i].y*fixed[i].y)*probabilities.pt1[i];
	}

	result.sigma2 -= np*(mu_x.x*mu_x.x+mu_x.y*mu_x.y)+(b1[0]*A[0]+b1[1]*A[1]+b1[2]*A[2]+b1[3]*A[3]);
	result.sigma2 /= dim*np;
}

void CPDEstimator::MaximumStepAsNonrigidTransform(const std::vector<util::point2d>& fixed, const std::vector<util::point2d>& moving,
											   const Probabilities& probabilities, double sigma2, NonrigidParams& result){
    int dim = 2;
	int m = moving.size();
	double np = 0;
	
	double A[m*m], w[m*dim], b[m*dim];
	
	for(int i = 0; i < m; i++){
		double* Ai = A+m*i;
		for(int j = 0; j < m; j++){
			Ai[j] = result.m_g[i][j]*probabilities.p1[i];		//d(P1)*G
		}
	}
	
	for(int i = 0; i < m; i++){
		A[m*i+i] += result.m_lambda*sigma2;
	}
	
	for(int i = 0; i < m; i++){
		double* bi = b+dim*i;
		bi[0] = probabilities.px[i].x - probabilities.p1[i]*moving[i].x;
		bi[1] = probabilities.px[i].y - probabilities.p1[i]*moving[i].y;
	}
	
	matrix_invert_inplace(m, A);
	matrix_product(m, m, m, dim, A, b, w);
	
	for(int i = 0; i < m; i++){
		double* wi = w+dim*i;
		result.m_w[i][0] = wi[0];
		result.m_w[i][1] = wi[1];
	}
	
	for(size_t i = probabilities.pt1.size(); i--;){	
		np += probabilities.pt1[i];
	}
	
	result.sigma2 = 0;
	for(size_t i = fixed.size(); i--;){			//for(size_t i = 0; i < fixed.size(); i++){
		result.sigma2 += (fixed[i].x*fixed[i].x+fixed[i].y*fixed[i].y)*probabilities.pt1[i];
	}
	
	std::vector<util::point2d> drifted;
	UpdatePointSetDrift(moving, result, drifted);
	
	for(size_t i = drifted.size(); i--;){			//for(size_t i = 0; i < fixed.size(); i++){
		result.sigma2 += (drifted[i].x*drifted[i].x+drifted[i].y*drifted[i].y)*probabilities.p1[i];
	}
	
	for(size_t i = drifted.size(); i--;){			//for(size_t i = 0; i < fixed.size(); i++){
		result.sigma2 -= 2*(drifted[i].x*probabilities.px[i].x+drifted[i].y*probabilities.px[i].y);
	}
	
	result.sigma2 /= dim*np;
}

void CPDEstimator::UpdatePointSet(const std::vector< util::point2d >& ori, const CPDEstimator::Params& params, std::vector< util::point2d >& result)
{
	result.resize(ori.size());
	const double* R = &(params.p[0]);
	const double* t = &(params.p[4]);
	
	for(size_t i = ori.size(); i--;){			//for(size_t i = 0; i < ori.size(); i++){
		result[i].x = (ori[i].x*R[0]+ori[i].y*R[1])+t[0];//params.scale*(ori[i].x*R[0]+ori[i].y*R[1])+t[0];
		result[i].y = (ori[i].x*R[2]+ori[i].y*R[3])+t[1];//params.scale*(ori[i].x*R[2]+ori[i].y*R[3])+t[1];
	}
}

void CPDEstimator::UpdatePointSet(const std::vector< util::point2d >& ori, const cv::Mat& A, std::vector< util::point2d >& result)
{
	result.resize(ori.size());
	const double* xf = (double*)(A.data);
	
	for(size_t i = ori.size(); i--;){			//for(size_t i = 0; i < ori.size(); i++){
		result[i].x = (ori[i].x*xf[0]+ori[i].y*xf[1])+xf[2];//params.scale*(ori[i].x*xf[0]+ori[i].y*xf[1])+t[0];
		result[i].y = (ori[i].x*xf[3]+ori[i].y*xf[4])+xf[5];//params.scale*(ori[i].x*xf[2]+ori[i].y*xf[3])+t[1];
	}
}

void CPDEstimator::FromCorrespondenceToMatch(const std::vector< util::point2d >& fixed, const std::vector< util::point2d >& moving, 
											 std::vector< int >& correspondence, std::vector< std::pair< int, int > >& rawmatches, double dist_err_tol) const
{
	float dist_err_thre = dist_err_tol*dist_err_tol;
	rawmatches.clear();
	for(int i = 0; i < correspondence.size(); i++){
		if(PT2DISTANCE(fixed[correspondence[i]], moving[i]) < dist_err_thre){
			rawmatches.push_back(std::make_pair(correspondence[i], i));
		}
	}
}

double CPDEstimator::GetMaxDiscernibleSigma(const std::vector<util::point2d>& pts1, const std::vector<util::point2d>& pts2)
{	
	SearchTree ref(pts1);
	int* idices = ref.IndicesArray();
	double* dists = ref.DistanceArray();
	double avg_dist = 0;
	
	for(int i = pts2.size(); i--;){
		ref.GetKNearestPoints(pts2[i], 2, idices, dists);
		if(&pts1 == &pts2){
			avg_dist += sqrt(dists[1]);
		}
		else{
			avg_dist += sqrt(dists[0]);
		}
	}
	avg_dist /= pts2.size();
	
	return avg_dist*.333333333333;//.2;		//2.5*sigma*2 = avg_dist or 1.5*sigma*2 = avg_dist
}

inline static void UpdatePointSetByCoarseParam(const std::vector<util::point2d>& moving, double alpha, double t0, double t1, std::vector<util::point2d>& nmoving)
{
	float calpha = cos(alpha);
	float salpha = sin(alpha);
	float A[4];
	A[0] = calpha; A[1] = salpha; A[2] = -salpha; A[3] = calpha;
	util::point2d tpt;
	
	for(int idx = nmoving.size(); idx--;){
		tpt.x = moving[idx].x+t0;
		tpt.y = moving[idx].y+t1;
		nmoving[idx].x = tpt.x*A[0]+tpt.y*A[1];
		nmoving[idx].y = tpt.x*A[2]+tpt.y*A[3];
	}
}

inline void CPDEstimator::AffineParamsProduct(double h, double alpi1, double alpi2, double t0, double t1, Params& tparams)
{
	double hmx[] = {h, 0, 0, 1};
	double alpi1mx[4], alpi2mx[4];
	alpi1mx[3] = alpi1mx[0] = cos(alpi1);
	alpi1mx[2] = -(alpi1mx[1] = sin(alpi1)); 
	alpi2mx[3] = alpi2mx[0] = cos(alpi2);
	alpi2mx[2] = -(alpi2mx[1] = sin(alpi2));
	double tmp[4];
	matrix_product(2, 2, 2, 2, hmx, alpi2mx, tmp);
	matrix_product(2, 2, 2, 2, alpi1mx, tmp, tparams.p);
	tparams.p[4] = tparams.p[0]*t0+tparams.p[1]*t1;
	tparams.p[5] = tparams.p[2]*t0+tparams.p[3]*t1;
}

void CPDEstimator::InitializeNonrigidParams()
{
	nr_params.m_beta = DEFAULT_BETA;
	nr_params.m_lambda = DEFAULT_LAMBDA;
	
	double k = -0.5/(nr_params.m_beta*nr_params.m_beta);
	
	nr_params.m_g.resize(moving.size());
	for(int i = nr_params.m_g.size(); i--;){
		nr_params.m_g[i].resize(moving.size());
		for(int j = 0; j < nr_params.m_g[i].size(); j++){
			nr_params.m_g[i][j] = exp(k*util::L2(moving[i], moving[j]));
		}
	}
	
	nr_params.m_w.resize(moving.size());
	for(int i = nr_params.m_w.size(); i--;){
		nr_params.m_w[i].resize(2);
		nr_params.m_w[i][0] = 0;
		nr_params.m_w[i][1] = 0;
	}
}

void CPDEstimator::ModifyProbabilities(const NonrigidParams& nr_params, Probabilities& probabilities)
{
	std::vector<double> wsum0(nr_params.m_w.size()), wsum1(nr_params.m_w.size());
	
	for(int j = 0; j < nr_params.m_g.size(); j++){			//gij
		wsum0[j] = 0; wsum1[j] = 0;
		for(int i = 0; i < nr_params.m_w.size(); i++){
			wsum0[j] += nr_params.m_w[i][0]*nr_params.m_g[i][j];
			wsum1[j] += nr_params.m_w[i][1]*nr_params.m_g[i][j];
		}
	}
	
	double wgw11 = 0, wgw22 = 0;		//mamtlab index, i.e., 0 -> 0+1
	for(int i = 0; i < nr_params.m_w.size(); i++){
		wgw11 += wsum0[i]*nr_params.m_w[i][0];
		wgw22 += wsum1[i]*nr_params.m_w[i][1];
	}
	
    probabilities.l += nr_params.m_lambda/2.0*(wgw11+wgw22);
}

void CPDEstimator::UpdatePointSetDrift(const std::vector< util::point2d >& ori, const NonrigidParams& params, std::vector< util::point2d >& result)
{
	result = ori;
	for(size_t i = ori.size(); i--;){
		for(int j = 0; j < params.m_g.size(); j++){
			result[i].x += params.m_g[i][j]*params.m_w[j][0];
			result[i].y += params.m_g[i][j]*params.m_w[j][1];
		}
	}
}		//result.points = moving + m_g * w;

bool CPDEstimator::GaussianMixedModelMatch(double dist_err_tol, std::vector< std::pair< int, int > >& matches, bool eigenlimit)
{
	double sigf = GetMaxDiscernibleSigma(fixed, fixed);
	double sigm = GetMaxDiscernibleSigma(moving, moving);
	default_sigma2 = pow(min(sigf, sigm), 2);
	
	InitializeNonrigidParams();
	
	Probabilities probabilities;
	points = moving;
	
	if (default_sigma2 == -1) {
		nr_params.sigma2 = DefaultSigma2(fixed, moving);
	}
	else{
		nr_params.sigma2 = default_sigma2;
	}
	
	size_t iter = max_iter;
	double ntol = cpd_err_tol + 10.0;
	double l = 0.;
	
	while(iter && ntol > cpd_err_tol && nr_params.sigma2 > 10*__DBL_EPSILON__){
		ComputeDirectGaussTransform(fixed, points, nr_params.sigma2, outliers_weight, probabilities);
		
		ntol = abs((probabilities.l - l)/probabilities.l);
		l = probabilities.l;

		MaximumStepAsNonrigidTransform(fixed, moving, probabilities, nr_params.sigma2, nr_params);
// 				MaximumStepAsRigidTransform(fixed, moving, probabilities, params.sigma2, params);
		UpdatePointSetDrift(moving, nr_params, points);
		iter--;
	}
	
// 	PrintPointSet(points);
	
	CalculateCorrespondence(fixed, points, fin_params.sigma2, outliers_weight, correspondence);
	
	if(eigenlimit){	
		double evec[4], eval[2];
		dgeev_driver(2, fin_params.p, evec, eval);
		if(eval[0] < 0.8 || eval[1] < 0.8 || eval[0] > 1.25 || eval[1] > 1.25){		// it is not a correct transform
// 			step_ratio *= 0.5;
			return false;
		}
	}
	
	FromCorrespondenceToMatch(fixed, points, correspondence, matches, dist_err_tol);		//a relax of threshold
		
	// 	std::cout<<matches.size()<<std::endl;	
	return true;
}

void CPDEstimator::MiniOptimizeByAffineGMM(const std::vector< util::point2d >& fixed, std::vector< util::point2d >& moving)
{
	Probabilities probabilities;
	points = moving;
	size_t iter = max_iter;
	double ntol = cpd_err_tol + 10.0;
	double l = 0.;
	
	fin_params.sigma2 = DefaultSigma2(fixed, moving);
	
// 	PrintPointSet(fixed);
// 	PrintPointSet(moving);
	
	while(iter && ntol > cpd_err_tol && fin_params.sigma2 > 10*__DBL_EPSILON__){
		ComputeDirectGaussTransform(fixed, points, fin_params.sigma2, outliers_weight, probabilities);
		
		ntol = abs((probabilities.l - l)/probabilities.l);
		l = probabilities.l;

		MaximumStepAsAffineTransform(fixed, moving, probabilities, fin_params.sigma2, fin_params);
// 			MaximumStepAsRigidTransform(fixed, moving, probabilities, params.sigma2, params);
		UpdatePointSet(moving, fin_params, points);
		iter--;
	}
	
	if(isnan(fin_params.sigma2)){
		return;
	}
	
	moving = points;
}

void CPDEstimator::MiniOptimizeByGMM(const std::vector< util::point2d >& fixed, const std::vector< util::point2d >& moving,
									 double dist_err_tol, std::vector< std::pair< int, int > >& matches, bool eigenlimit)
{
	Probabilities probabilities;
	points = moving;
	size_t iter = max_iter;
	double ntol = cpd_err_tol + 10.0;
	double l = 0.;

	double sigf = GetMaxDiscernibleSigma(fixed, fixed);
	double sigm = GetMaxDiscernibleSigma(moving, moving);
	fin_params.sigma2 = pow(min(sigf, sigm), 2);
	matches.clear();
// 	fin_params.sigma2 = default_sigma2;

// 	PrintPointSet(fixed);
// 	PrintPointSet(moving);

	while(iter && ntol > cpd_err_tol && fin_params.sigma2 > 10*__DBL_EPSILON__){
		ComputeDirectGaussTransform(fixed, points, fin_params.sigma2, outliers_weight, probabilities);

		ntol = abs((probabilities.l - l)/probabilities.l);
		l = probabilities.l;

		MaximumStepAsAffineTransform(fixed, moving, probabilities, fin_params.sigma2, fin_params);
// 			MaximumStepAsRigidTransform(fixed, moving, probabilities, params.sigma2, params);
		UpdatePointSet(moving, fin_params, points);
		iter--;
	}

	if(isnan(fin_params.sigma2)){
		return;
	}

// 	PrintPointSet(points);

	CalculateCorrespondence(fixed, points, fin_params.sigma2, outliers_weight, correspondence);

	double evec[4], eval[2];
	dgeev_driver(2, fin_params.p, evec, eval);

	if(eigenlimit){
		if(eval[0] < 0.8 || eval[1] < 0.8 || eval[0] > 1.25 || eval[1] > 1.25){		// it is not a correct transform
			return;
		}
	}

	FromCorrespondenceToMatch(fixed, points, correspondence, matches, dist_err_tol);
}

// CPDEstimator::CPDEstimator(const std::vector<util::point2d>& __fids1, const std::vector<util::point2d>& __fids2): ori_fixed(__fids1), ori_moving(__fids2),
// 						fixed(__fids1), moving(__fids2), default_sigma2(DEFAULT_SIGMA2), cpd_err_tol(ERROR_TOL), outliers_weight(OUTLIERS_WEIGHT), max_iter(MAX_CPD_ITERATION)
// {
//
// }

CPDEstimator::CPDEstimator(const std::vector<util::point2d>& __fids1, const std::vector<util::point2d>& __fids2, int method): ori_fixed(__fids1), ori_moving(__fids2),
						fixed(__fids1), moving(__fids2), default_sigma2(DEFAULT_SIGMA2), cpd_err_tol(ERROR_TOL), outliers_weight(OUTLIERS_WEIGHT), max_iter(MAX_CPD_ITERATION)
{
	if(method)
	{
		Normalizer::GlobalNormalize(fixed, moving, fixed_center, moving_center, &scale);
		new (&ftree) SearchTree(fixed);
	}
}


void CPDEstimator::CalculateMatch(const std::vector< util::point2d >& fids1, const std::vector< util::point2d >& fids2,
									   std::vector< std::pair<int, int> >& rawmatches, std::vector<double>& dists, double dist_err_tol)
{
	SearchTree ref(fids1);
	int* idices = ref.IndicesArray();
	double* distances = ref.DistanceArray();
	float dist_err_thre = dist_err_tol*dist_err_tol;
	
	rawmatches.clear();
	dists.clear();
	
	for(int i = fids2.size(); i--;){
		ref.GetKNearestPoints(fids2[i], 1, idices, distances);
		if(distances[0] < dist_err_thre){
			rawmatches.push_back(std::make_pair(idices[0], i));
			dists.push_back(distances[0]);
		}
	}
}

/**
 * |a b|  |t0|
 * |c d|  |t1|

/**xf=params:[a b c d t0 t1]*/
void CPDEstimator::ComputeTransformByLeastSquares(const std::vector< util::point2d >& fixed, const std::vector< util::point2d >& moving, 
												  const std::vector< std::pair< int, int > >& matches, CPDEstimator::Params& params)
{
	size_t n = matches.size();
	int num_eqs = 2*n;
    int num_vars = 6;

    double *As = new double[num_eqs*num_vars];
    double *bs = new double[num_eqs];

    for(size_t i = n; i--;){
        double* A = As+12*i;
        double* b = bs+2*i;
		
        A[0] = moving[matches[i].second].x; A[1] = moving[matches[i].second].y; A[2] = 0; A[3] = 0; A[4] = 1; A[5] = 0;
		A[6] = 0; A[7] = 0; A[8] = moving[matches[i].second].x; A[9] = moving[matches[i].second].y; A[10] = 0; A[11] = 1;
        b[0] = fixed[matches[i].first].x;
        b[1] = fixed[matches[i].first].y;
    }
    
    dgelsy_driver(As, bs, params.p, num_eqs, num_vars, 1);

	delete [] As;
	delete [] bs;
}

static bool pairCompare(const std::pair<double, int>& firstElem, const std::pair<double, int>& secondElem){
  return firstElem.first < secondElem.first;
}

void CPDEstimator::CoarseAlignmentByGMM(const std::vector<util::point2d>& fixed, std::vector<util::point2d>& moving, double sigma2, CPDEstimator::Params& param, double step_ratio)
{
#define PI_180				0.0174532925199433
#define SPAN				0.65
	double grid_step = step_ratio*sqrt(sigma2);//2.5*sqrt(sigma2);		//grid_step = 2.5*sigma
	double angle_step = atan(grid_step);
	double corr, max_corr = -1;

	std::vector<util::point2d> nmoving(moving.size());
// 	EX_TRACE("coarse search grid step = %0.2lfpx, angle step = %0.2lf\n", grid_step*scale, angle_step/PI_180);

	for(double t0 = -SPAN; t0 < SPAN; t0 += grid_step){
		for(double t1 = -SPAN; t1 < SPAN; t1 += grid_step){
			for(double alpi = -0.523598775598299; alpi < 0.523598775598299; alpi += angle_step){		//30*pi/180 = 0.523598775598299
				UpdatePointSetByCoarseParam(moving, alpi, t0, t1, nmoving);
// 				corr = ComputeDirectGaussCorrelation(fixed, nmoving, sigma2);
				corr = ComputeTreeGaussCorrelation(fixed, nmoving, sigma2);
				if(corr > max_corr){
					max_corr = corr;
					param.p[3] = alpi;
					param.p[4] = t0;
					param.p[5] = t1;
				}
			}
		}
	}
// 	std::cout<<max_corr<<std::endl;
	UpdatePointSetByCoarseParam(moving, param.p[3], param.p[4], param.p[5], nmoving);
	moving = nmoving;
	double tmp[2];
	param.p[0] = cos(param.p[3]); param.p[1] = sin(param.p[3]);
	param.p[2] = -param.p[1]; param.p[3] = param.p[0];
	tmp[0] = param.p[0]*param.p[4]+param.p[1]*param.p[5];
	tmp[1] = param.p[2]*param.p[4]+param.p[3]*param.p[5];
	param.p[4] = tmp[0]; param.p[5] = tmp[1];
// 	PrintPointSet(moving);
}



bool CPDEstimator::GaussianMixedModelMatch_Global(double dist_err_tol, std::vector<std::pair<int, int>>& matches, bool eigenlimit)
{
	if((int)fixed.size() < MIN_PT_NUM || (int)fixed.size() < MIN_PT_NUM){
		EX_TRACE("Insufficient number of fiducial markers\n");
		return false;
	}

	double sigf = GetMaxDiscernibleSigma(fixed, fixed);
	double sigm = GetMaxDiscernibleSigma(moving, moving);
	default_sigma2 = pow(min(sigf, sigm), 2);
	double step_ratio = 2.5;

	bool solved = false;
	int count = 0;
	while(!solved && count < 2){
// 		PrintPointSet(fixed);
// 		PrintPointSet(moving);

		CoarseAlignmentByGMM(fixed, moving, default_sigma2, cos_params, step_ratio);

// 		PrintPointSet(moving);

		Probabilities probabilities;
		points = moving;

		if (default_sigma2 == -1) {
			fin_params.sigma2 = DefaultSigma2(fixed, moving);
		}
		else{
			fin_params.sigma2 = default_sigma2;
		}

		size_t iter = max_iter;
		double ntol = cpd_err_tol + 10.0;
		double l = 0.;

		while(iter && ntol > cpd_err_tol && fin_params.sigma2 > 10*__DBL_EPSILON__){
			ComputeDirectGaussTransform(fixed, points, fin_params.sigma2, outliers_weight, probabilities);

			ntol = abs((probabilities.l - l)/probabilities.l);
			l = probabilities.l;

			MaximumStepAsAffineTransform(fixed, moving, probabilities, fin_params.sigma2, fin_params);
// 				MaximumStepAsRigidTransform(fixed, moving, probabilities, params.sigma2, params);
			UpdatePointSet(moving, fin_params, points);
			iter--;
		}

		// 		PrintPointSet(points);

		CalculateCorrespondence(fixed, points, fin_params.sigma2, outliers_weight, correspondence);

		if(eigenlimit){
			double evec[4], eval[2];
			dgeev_driver(2, fin_params.p, evec, eval);
			if(eval[0] < 0.8 || eval[1] < 0.8 || eval[0] > 1.25 || eval[1] > 1.25){		// it is not a correct transform
				step_ratio *= 0.5;
				count++;
				continue;
			}
		}

		FromCorrespondenceToMatch(fixed, points, correspondence, matches, dist_err_tol*1.25);		//a relax of threshold
		if(matches.size() < 5 || (matches.size() < 9 && matches.size() < fixed.size()*COVERAGE_REFRESH_RATIO*.15 &&  matches.size() < moving.size()*COVERAGE_REFRESH_RATIO*.15)){
			step_ratio *= 0.5;		// in the case of incorrect start point
			if(step_ratio < 0.5){
				return false;		//we failed
			}
		}
		else{
			solved = true;
		}
	}
	// 	std::cout<<matches.size()<<std::endl;
	return true;
}


void CPDEstimator::GlobalMatch(double dist_err_tol, std::vector<std::pair<util::point2d, util::point2d>>& matchvec, std::vector<cv::Mat>& hset, bool test, bool eigenlimit)
{
	pts_dist_tol = dist_err_tol/scale;

	std::vector<std::pair<int, int> > matches;
	if(!GaussianMixedModelMatch_Global(pts_dist_tol, matches, eigenlimit)){
		return;
	}
// std::cout<<"???????????"<<std::endl;
	std::vector<util::point2d> nfixed = fixed;				//used for calculation
	std::vector<util::point2d> nmoving = moving;			//used for calculation
	std::vector<util::point2d> norif = ori_fixed;			//used for datasaving
	std::vector<util::point2d> norim = ori_moving;			//used for datasaving

	std::vector<std::vector< std::pair< util::point2d, util::point2d > > > pairsets;
	std::vector<Params> pvec;

	std::vector<std::pair< util::point2d, util::point2d > > prepairs;

	bool split_pointset = true; bool first_itr = true;
	double dist_err_thre = pts_dist_tol*pts_dist_tol;		//square distance error
	while(split_pointset){
		Params nparams;
		std::vector<std::pair<int, int> > indices;
		std::vector<std::pair<int, int> > premidices;
		std::vector<std::pair<int, int> > nmatches;

		if(matches.size() < 6){					//there are no matches any more
            if(!pairsets.size()){
				break;
			}
			pairsets[pairsets.size()-1].insert(pairsets[pairsets.size()-1].end(), prepairs.begin(), prepairs.end());
			break;
		}
		if(!first_itr){
			std::vector<util::point2d> tfixed(matches.size()), tmoving(matches.size());
			std::vector<int> tforidx(matches.size()), tmoridx(matches.size());
			for(int i = 0; i < matches.size(); i++){
				tfixed[i] = nfixed[matches[i].first];
				tforidx[i] = matches[i].first;
				tmoving[i] = nmoving[matches[i].second];
				tmoridx[i] = matches[i].second;
			}
			MiniOptimizeByGMM(tfixed, tmoving, pts_dist_tol, matches, eigenlimit);		//used to solve uncertain assign problem
			for(int i = 0; i < matches.size(); i++){
				matches[i].first = tforidx[matches[i].first];
				matches[i].second = tmoridx[matches[i].second];
			}
			if(matches.size() < 6){				//it indicates that the remains are with bad distribution
				break;
			}
		}
		ComputeTransformByLeastSquares(nfixed, nmoving, matches, nparams);	//!WORNING use GMMaffine?

		UpdatePointSet(nmoving, nparams, points);

		std::vector<double> dists;
		CalculateMatch(nfixed, points, matches, dists, pts_dist_tol*2);	//get all the possible matches

// 		std::cout<<matches.size()<<std::endl;

		std::vector<std::pair<double, int> > candidates;
		for(int i = 0; i < matches.size(); i++){
			candidates.push_back(std::make_pair(dists[i], i));
		}

		std::sort(candidates.begin(), candidates.end(), pairCompare);

		if(!candidates.size()){
			break;
		}
		else if(candidates[candidates.size()-1].first < dist_err_thre){		//the two point sets have no warp
			indices = matches;
			pvec.push_back(nparams);
			split_pointset = false;
			prepairs.clear();
		}
		else{
			prepairs.clear();

			int candix;
			double warp_dist_err_thre = pts_dist_tol*pts_dist_tol*.75*.75;

			for(candix = candidates.size(); candix--;){
				if(candidates[candix].first < warp_dist_err_thre){			//find the point pairs that are close enough
					break;
				}
				if(candidates[candix].first < dist_err_thre){				//save the matches between thre*.75 and thre in the case that next routine failed
					premidices.push_back(matches[candidates[candix].second]);
				}
			}
			if(candix < 3){
				break;
			}

// 			std::cout<<"premidices: "<<premidices.size()<<std::endl;
// 			std::cout<<"candix: "<<candix<<std::endl;

			for(int i = candix+1; i < candidates.size(); i++){				//next turn's seed;
				nmatches.push_back(matches[candidates[i].second]);
			}

			for(int idx = candix+1; idx--;){
				indices.push_back(matches[candidates[idx].second]);
			}
			pvec.push_back(nparams);
		}

		std::vector< std::pair< util::point2d, util::point2d > > pairs;
		for(int i = 0; i < indices.size(); i++){
			pairs.push_back(std::make_pair(norif[indices[i].first], norim[indices[i].second]));
		}

		pairsets.push_back(pairs);

		for(int i = 0; i < premidices.size(); i++){
			prepairs.push_back(std::make_pair(norif[premidices[i].first], norim[premidices[i].second]));
		}

		std::vector<util::point2d> tfixed, torif;
		for(int i = 0; i < nfixed.size(); i++){
			bool find = false;
			for(int j = 0; j < indices.size(); j++){
				if(indices[j].first == i){
					find = true;
					break;
				}
			}
			if(!find){
				tfixed.push_back(nfixed[i]);
				torif.push_back(norif[i]);
			}
			for(int j = 0; j < nmatches.size(); j++){
				if(nmatches[j].first == i){
					nmatches[j].first = tfixed.size()-1;
				}
			}
		}
		nfixed = tfixed; norif = torif;

		std::vector<util::point2d> tmoving, torim;
		for(int i = 0; i < nmoving.size(); i++){
			bool find = false;
			for(int j = 0; j < indices.size(); j++){
				if(indices[j].second == i){
					find = true;
					break;
				}
			}
			if(!find){
				tmoving.push_back(nmoving[i]);
				torim.push_back(norim[i]);
			}
			for(int j = 0; j < nmatches.size(); j++){
				if(nmatches[j].second == i){
					nmatches[j].second = tmoving.size()-1;
				}
			}
		}
		nmoving = tmoving; norim = torim;

		matches = nmatches;
		first_itr = false;

	}

	matchvec.clear();
	for(int i = 0; i < pairsets.size(); i++){
		std::vector<std::pair<util::point2d, util::point2d> >& mm = pairsets[i];
		matchvec.insert(matchvec.end(), mm.begin(), mm.end());
	}

	hset.clear();
	double A[4], t[2];
	for(int i = 0; i < pvec.size(); i++){			// change the transform matrix according to the coordinate of ori point sets
		A[0] = pvec[i].p[0]*cos_params.p[0]+pvec[i].p[1]*cos_params.p[2];
		A[1] = pvec[i].p[0]*cos_params.p[1]+pvec[i].p[1]*cos_params.p[3];
		A[2] = pvec[i].p[2]*cos_params.p[0]+pvec[i].p[3]*cos_params.p[2];
		A[3] = pvec[i].p[2]*cos_params.p[1]+pvec[i].p[3]*cos_params.p[3];
		t[0] = pvec[i].p[4];
		t[1] = pvec[i].p[5];
		t[0] += pvec[i].p[0]*cos_params.p[4]+pvec[i].p[1]*cos_params.p[5];
		t[1] += pvec[i].p[2]*cos_params.p[4]+pvec[i].p[3]*cos_params.p[5];
		t[0] *= scale; t[1] *= scale;
		t[0] += -A[0]*moving_center.x-A[1]*moving_center.y+fixed_center.x;
		t[1] += -A[2]*moving_center.x-A[3]*moving_center.y+fixed_center.y;
		memcpy(pvec[i].p, A, sizeof(double)*4);
		memcpy(pvec[i].p+4, t, sizeof(double)*2);
// std::cout<<"A"<<A[0]<<"    "<<A[1]<<std::endl;
// 		CvMat* h = cvCreateMat(3, 3, CV_64FC1);
		cv::Mat h =cv::Mat::zeros(3, 3, CV_64FC1);
		h.at<double>(0,0)=A[0];h.at<double>(0,1)=A[1];h.at<double>(0,2)=t[0];
		h.at<double>(1,0)=A[2];h.at<double>(1,1)=A[3];h.at<double>(1,2)=t[1];
		h.at<double>(2,0)=0;h.at<double>(2,1)=0;h.at<double>(2,2)=1;
// 		cvmSet(h,0,0,A[0]); cvmSet(h,0,1,A[1]); cvmSet(h,0,2,t[0]);
// 		cvmSet(h,1,0,A[2]); cvmSet(h,1,1,A[3]); cvmSet(h,1,2,t[1]);
// 		cvmSet(h,2,0,0); cvmSet(h,2,1,0); cvmSet(h,2,2,1);
// std::cout<<"h:"<<h.at<double>(0,1)<<"   "<<h.at<double>(1,0)<<"   "<<A[1]<<"   "<<A[2]<<std::endl;
		cv::invert(h, h);			//compatible with random 4 points methods

		hset.push_back(h);
// 		std::cout<<"h:"<<h<<std::endl;
	}


}


void CPDEstimator::PointDriftEstimation(double dist_err_tol, std::vector< std::pair< util::point2d, util::point2d > >& matchvec, 
										std::vector<cv::Mat >& hset, bool test, bool eigenlimit)
{
	if(fixed.size() < 2){
		EX_TRACE("Insufficient number of fiducial markers\n");
		return;
	}
	
	if(hset.size()){			//correct moving according to affine transformation
		std::vector<util::point2d> moving_trans;
		initial_trans_mat = hset[0];
		UpdatePointSet(moving, initial_trans_mat, moving_trans);
// 		hset[0].Print(std::cout);
		moving = moving_trans;
	}
	
	Normalizer::GlobalNormalize(fixed, moving, fixed_center, moving_center, &scale);
	pts_dist_tol = dist_err_tol/scale;
		
// 	PrintPointSet(fixed);
// 	PrintPointSet(moving);
	
	MiniOptimizeByAffineGMM(fixed, moving);
	
	new (&ftree) SearchTree(fixed);
	
	std::vector<std::pair<int, int> > matches;
	
	GaussianMixedModelMatch(.85*pts_dist_tol, matches, eigenlimit);
	
	for(int i = matches.size(); i--;){
// 		matchvec.push_back(std::make_pair(ori_fixed[matches[i].first], ori_moving[matches[i].second]));
		matchvec.push_back(std::make_pair(ori_moving[matches[i].second], ori_fixed[matches[i].first]));	//keep the order as 4p method
	}

}


void CPDEstimator::DistributionCorrection(double dist_err_tol, std::vector< std::pair< util::point2d, util::point2d > >& matchvec)
{
	SearchTree ref1(ori_fixed);
	SearchTree ref2(ori_moving);
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
