
	template <typename T>
ControllerADP<T>::ControllerADP(const SymmetricMatrix& Q, const SymmetricMatrix& R, const double delta,Step* stepf)
	:mn(Q.size()[0]), mm(R.size()[0]), mK0(mn,mm), mKadp(mn,mm), mQ(&Q), mR(&R), mP(mn), 
	mdelta(delta), itx(mxx.begin()), itxx(mIxx.begin()), itxu(mIxu.begin()),
	mThetaInv(mm*mn+mn*(mn+1)/2, mm*mn+mn*(mn+1)/2), 
	mBigV(mm*mn+mn*(mn+1)/2,0),
	mBigr(new Matrix(mn*mn,mm*mn+mn*(mn+1)/2)),
	mADPalg(nullptr),
	mResult(nullptr)
{
	std::shared_ptr<AlgorithmADP> tem (new T());
	mADPalg = (tem->Creat(*mQ,*mR,mP,mK0,stepf));
	//long double eps = 1e-10;
	//mThetaInv=mThetaInv+1/eps; 
	mThetaInv=mThetaInv+1e10;
}
	template <typename T>
ControllerADP<T>::ControllerADP(const SymmetricMatrix& Q, const SymmetricMatrix& R, const double delta, const SymmetricMatrix& P, Step* stepf)
	:mn(Q.size()[0]), mm(R.size()[0]), mK0(mn,mm), mKadp(mn,mm), mQ(&Q), mR(&R), mP(P), 
	mdelta(delta), itx(mxx.begin()), itxx(mIxx.begin()), itxu(mIxu.begin()),
	mThetaInv(mm*mn+mn*(mn+1)/2, mm*mn+mn*(mn+1)/2), 
	mBigV(mm*mn+mn*(mn+1)/2,0),
	mBigr(new Matrix(mn*mn,mm*mn+mn*(mn+1)/2)),
	mADPalg(nullptr),
	mResult(nullptr)
{
	std::shared_ptr<AlgorithmADP> tem (new T());
	mADPalg = (tem->Creat(*mQ,*mR,mP,mK0,stepf));
	//long double eps = 1e-10;
	//mThetaInv=mThetaInv+1/eps; 
	mThetaInv=mThetaInv+1e10;
}
	template <typename T>
ControllerADP<T>::ControllerADP(const SymmetricMatrix& Q, const SymmetricMatrix& R, const double delta, const Matrix& K, Step* stepf)
	:mn(Q.size()[0]), mm(R.size()[0]), mK0(K), mKadp(K), mQ(&Q), mR(&R), mP(mn), 
	mdelta(delta), itx(mxx.begin()), itxx(mIxx.begin()), itxu(mIxu.begin()),
	mThetaInv(mm*mn+mn*(mn+1)/2, mm*mn+mn*(mn+1)/2), 
	mBigV(mm*mn+mn*(mn+1)/2,0),
	mBigr(new Matrix(mn*mn,mm*mn+mn*(mn+1)/2)),
	mADPalg(nullptr),
	mResult(nullptr)
{

	std::shared_ptr<AlgorithmADP> tem (new T());
	mADPalg = (tem->Creat(*mQ,*mR,mP,mK0,stepf));
	//long double eps = 1e-10;
	//mThetaInv=mThetaInv+1/eps; 
	mThetaInv=mThetaInv+1e10;
}

	template <typename T>
ControllerADP<T>::ControllerADP(const SymmetricMatrix& Q, const SymmetricMatrix& R, const double delta, const SymmetricMatrix& P, const Matrix& K, Step* stepf)
	:mn(Q.size()[0]), mm(R.size()[0]), mK0(K), mKadp(K), mQ(&Q), mR(&R), mP(P), mdelta(delta), itx(mxx.begin()), itxx(mIxx.begin()), itxu(mIxu.begin()),
	mThetaInv(mm*mn+mn*(mn+1)/2, mm*mn+mn*(mn+1)/2), 
	mBigV(mm*mn+mn*(mn+1)/2,0),
	mBigr(new Matrix(mn*mn,mm*mn+mn*(mn+1)/2)),
	mADPalg(nullptr),
	mResult(nullptr)
{
	std::shared_ptr<AlgorithmADP> tem (new T());
	mADPalg = (tem->Creat(*mQ,*mR,mP,mK0,stepf));
	//long double eps = 1e-10;
	mThetaInv=mThetaInv+1e10;
	//mThetaInv=mThetaInv+1/eps; 
}

	template <>
const Matrix& ControllerADP<AlgorithmPI>::learner(const std::vector<double>& x, const std::vector<double>& u, const double dt, const double t)
{
	if(mxx.size()==0) 
	{
		mxx.push_back(vecs(prod(x,x)));
		++itx;
	}else mxx.push_back(vecs(prod(x,x)));


	if(mIxx.size()==0) 
	{
		mIxx.push_back(-dt*vec(kProd(x,x)));
		++itxx;
	}
	else mIxx.push_back(mIxx.back()-dt*vec(kProd(x,x)));




	if(mIxu.size()==0) 
	{
		mIxu.push_back(-2*dt*vec(kProd(x,*mR*u)));
		++itxu;
	} else mIxu.push_back(mIxu.back()-2*dt*vec(kProd(x,*mR*u)));

	if (t > mdelta && mxx.size()>1 && mIxx.size()>1 && mIxu.size()>1){
		//std::cout << mdelta << std::endl;
		mdelta +=500*dt;
		itx = mxx.erase(itx);
		itxx = mIxx.erase(itxx);
		itxu = mIxu.erase(itxu);

		// old one
		//std::vector<double> vec0 = mIxx.back() - *itxx;
		////std::vector<double> vec1 = mIxx.back() - *itxx;
		//std::vector<double> vec1;
		//vec1.reserve(itxx->size()+itx->size());
		//vec1 = mxx.back() - *itx;
		//std::vector<double> vec2 = vec(2*kProd(Diagonal(mQ.size()[0]),mR*mKadp)*vec0) + mIxu.back() - *itxu;
		//vec1.insert( vec1.end(), vec2.begin(), vec2.end() );
		////disp(vec1);
		//LS(vec1, vec0*vec(mQ+mKadp.t()*mR*mKadp));


		// new one
		//std::vector<double> vec0 = (mxx.back() - *itx)/(x*x);
		//std::vector<double> vec1 = (mIxx.back() - *itxx)/(x*x);
		//std::vector<double> vec2 = (mIxu.back() - *itxu)/(x*x);

		mResult = &mADPalg->online((mxx.back() - *itx)/(x*x), (mIxx.back() - *itxx)/(x*x), (mIxu.back() - *itxu)/(x*x), mBigr, mThetaInv, mBigV);
		//double err= optResult[0].F();
		mP = (*mResult)[1];
		mKadp= (*mResult)[2];
		//disp(mP);
		//mBigV = vec(mThetaInv * *mBigr * vec(mP));

		//std::cout << "the error is " << err << std::endl;
		//mP.disp();


		itx = mxx.end();
		itxx = mIxx.end();
		itxu = mIxu.end();
		--itx;
		--itxx;
		--itxu;

	}

	return mKadp;
}


	template <typename T>
const Matrix& ControllerADP<T>::learner(const std::vector<double>& x, const std::vector<double>& u, const double dt, const double t)
{
	if(mxx.size()==0) 
	{
		mxx.push_back(vec(kProd(x,x)));
		++itx;
	}else mxx.push_back(vec(kProd(x,x)));

	if(mIxx.size()==0) 
	{
		mIxx.push_back(dt*vecs(prod(x,x)));
		++itxx;
	}
	else mIxx.push_back(mIxx.back()+dt*vecs(prod(x,x)));

	if(mIxu.size()==0) 
	{
		mIxu.push_back(2*dt*vec(kProd(x,*mR*u)));
		++itxu;
	} else mIxu.push_back(mIxu.back()+2*dt*vec(kProd(x,*mR*u)));



	if (t > mdelta && mxx.size()>1 && mIxx.size()>1 && mIxu.size()>1){
		//std::cout << mdelta << std::endl;
		mdelta += 500 * dt;
		itx = mxx.erase(itx);
		itxx = mIxx.erase(itxx);
		itxu = mIxu.erase(itxu);

		// this is the old one
		//std::vector<double> vec0 = mxx.back() - *itx;
		////std::vector<double> vec1 = mIxx.back() - *itxx;
		//std::vector<double> vec1;
		//vec1.reserve(itxx->size()+itxu->size());
		//vec1 = mIxx.back() - *itxx;
		//std::vector<double> vec2 = mIxu.back() - *itxu;
		//vec1.insert( vec1.end(), vec2.begin(), vec2.end() );
		//*mBigr= *mBigr + prod(vec1,vec0);// half online case;
		////std::cout<< double(vec0*vec(mP)) << std::endl;
		////disp(vec1);

		//LS(vec1, double(vec0*vec(mP)));

		////if((mBigV-mBigVold)*(mBigV-mBigVold)<1e-10)
		////{
		//// online method
		//std::vector<Matrix> optResult = mADPalg->onlineI(mBigV);
		//mP = optResult[1];
		//mKadp= optResult[2];
		//mBigV = vec(mThetaInv * *mBigr * vec(mP));

		//std::cout<< "|x| = " << norm(x) << std::endl;



		// this is the new one
		//std::vector<double> vec0 = (mxx.back() - *itx)/(x*x);
		//std::vector<double> vec1 = (mIxx.back() - *itxx)/(x*x);
		//std::vector<double> vec2 = (mIxu.back() - *itxu)/(x*x);

		//std::vector<Matrix> optResult = mADPalg->online(vec0, vec1, vec2, mBigr, mThetaInv, mBigV);
		//std::vector<Matrix> optResult = mADPalg->online((mxx.back() - *itx)/(x*x), (mIxx.back() - *itxx)/(x*x), (mIxu.back() - *itxu)/(x*x), mBigr, mThetaInv, mBigV);
		mResult = &mADPalg->online((mxx.back() - *itx)/(x*x), (mIxx.back() - *itxx)/(x*x), (mIxu.back() - *itxu)/(x*x), mBigr, mThetaInv, mBigV);
		//double err= norm(optResult[0]);
		mP = (*mResult)[1];
		mKadp= (*mResult)[2];
		//mBigV = vec(mThetaInv * *mBigr * vec(mP)); // half online

		//std::cout << "the error is " << err << std::endl;
		//mP.disp();

		if(norm((*mResult)[0])<1e-10)                                                                            // convergence
		{
			mThetaInv=mThetaInv*0+1e10;
			mBigV = mBigV*0;
			mxx.clear();
			mIxx.clear();
			mIxu.clear();
			itx = mxx.begin();
			itxx = mIxx.begin();
			itxu = mIxu.begin();
			//mBigr.reset(new Matrix(mn*mn,mm*mn+mn*(mn+1)/2));
			*mBigr = *mBigr * 0;
			mADPalg->resetStep();
		}else
		{
			itx = mxx.end();
			itxx = mIxx.end();
			itxu = mIxu.end();
			--itx;
			--itxx;
			--itxu;
		}



		// half online half offline method
		//for (int k=1;k<=2000;++k){
		//std::vector<Matrix> optResult = mADPalg->onlineI(mBigV, mQ, mR);
		//mP = optResult[1];
		//mKadp= optResult[2];
		//mP.disp();
		//mBigV = vec(mThetaInv * *mBigr * vec(mP));
		//}

		//}
	}
	return mKadp;
}


	template <typename T>
const std::vector<double> ControllerADP<T>::input(const std::vector<double>& x, const double dt,  const double t, noise noif)
{	
	double mag = 1/(norm(x)+10);
	std::vector<double> u = vec(-mKadp*x)+noif(mag,mm,t);
	learner(x,u,dt,t);

	return u;

}

	template <>
const std::vector<double> ControllerADP<AlgorithmPI>::input(const std::vector<double>& x, const double dt,  const double t, noise noif)
{	
	double mag = 1/(norm(x)+10);
	std::vector<double> u = vec(-mKadp*x)+noif(mag,mm,t);
	learner(x,u,dt,t);

	if (dt<1e-20) 
		return noif(mag,mm,t);

	return u;

}



template <typename T>
void ControllerADP<T>::dispAll() const 
{
	//std::cout <<"what is in mxx" << std::endl;
	//for (auto itx = mxx.begin();itx!=mxx.end();++itx)
	//disp(*itx);
	//std::cout <<"what is in mIxx" << std::endl;
	//for (auto itxx = mIxx.begin();itxx!=mIxx.end();++itxx)
	//disp(*itxx);
	//std::cout <<"what is in mIxu" << std::endl;
	//for (auto itxu = mIxu.begin();itxu!=mIxu.end();++itxu)
	//disp(*itxu);
	std::cout <<"what is Kadp" << std::endl;
	disp(mKadp);
	std::cout <<"what is K0" << std::endl;
	disp(mK0);
	std::cout <<"what is P" << std::endl;
	disp(mP);

}

