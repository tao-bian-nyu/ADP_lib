
	template <typename T>
	ControllerADP<T>::ControllerADP(const SymmetricMatrix& Q, const SymmetricMatrix& R, const double delta,Step* stepf)
:mK0(Q.size()[0],R.size()[0]), mKadp(Q.size()[0],R.size()[0]), mQ(Q), mR(R), mP(Q.size()[0]), 
	mdelta(delta), itx(mxx.begin()), itxx(mIxx.begin()), itxu(mIxu.begin()),
	mThetaInv(R.size()[0]*Q.size()[0]+Q.size()[0]*(Q.size()[0]+1)/2, R.size()[0]*Q.size()[0]+Q.size()[0]*(Q.size()[0]+1)/2), 
	mBigV(R.size()[0]*Q.size()[0]+Q.size()[0]*(Q.size()[0]+1)/2,0),
	mBigTheta(new Matrix(R.size()[0]*Q.size()[0]+Q.size()[0]*(Q.size()[0]+1)/2,R.size()[0]*Q.size()[0]+Q.size()[0]*(Q.size()[0]+1)/2)),
	mBigr(new Matrix(Q.size()[0]*Q.size()[0],R.size()[0]*Q.size()[0]+Q.size()[0]*(Q.size()[0]+1)/2)),
	mADPalg(nullptr)
{
	std::shared_ptr<AlgorithmADP> tem (new T());
	mADPalg = (tem->Creat(mQ,mR,mP,mK0,stepf));
	long double eps = 1e-10;
	mThetaInv=mThetaInv+1/eps; 
}
	template <typename T>
	ControllerADP<T>::ControllerADP(const SymmetricMatrix& Q, const SymmetricMatrix& R, const double delta, const SymmetricMatrix& P, Step* stepf)
:mK0(Q.size()[0],R.size()[0]), mKadp(Q.size()[0],R.size()[0]), mQ(Q), mR(R), mP(P), 
	mdelta(delta), itx(mxx.begin()), itxx(mIxx.begin()), itxu(mIxu.begin()),
	mThetaInv(R.size()[0]*Q.size()[0]+Q.size()[0]*(Q.size()[0]+1)/2, R.size()[0]*Q.size()[0]+Q.size()[0]*(Q.size()[0]+1)/2), 
	mBigV(R.size()[0]*Q.size()[0]+Q.size()[0]*(Q.size()[0]+1)/2,0),
	mBigTheta(new Matrix(R.size()[0]*Q.size()[0]+Q.size()[0]*(Q.size()[0]+1)/2,R.size()[0]*Q.size()[0]+Q.size()[0]*(Q.size()[0]+1)/2)),
	mBigr(new Matrix(Q.size()[0]*Q.size()[0],R.size()[0]*Q.size()[0]+Q.size()[0]*(Q.size()[0]+1)/2)),
	mADPalg(nullptr)
{
	std::shared_ptr<AlgorithmADP> tem (new T());
	mADPalg = (tem->Creat(mQ,mR,mP,mK0,stepf));
	long double eps = 1e-10;
	mThetaInv=mThetaInv+1/eps; 
}
	template <typename T>
	ControllerADP<T>::ControllerADP(const SymmetricMatrix& Q, const SymmetricMatrix& R, const double delta, const Matrix& K, Step* stepf)
:mK0(K), mKadp(K), mQ(Q), mR(R), mP(Q.size()[0]), 
	mdelta(delta), itx(mxx.begin()), itxx(mIxx.begin()), itxu(mIxu.begin()),
	mThetaInv(R.size()[0]*Q.size()[0]+Q.size()[0]*(Q.size()[0]+1)/2, R.size()[0]*Q.size()[0]+Q.size()[0]*(Q.size()[0]+1)/2), 
	mBigV(R.size()[0]*Q.size()[0]+Q.size()[0]*(Q.size()[0]+1)/2,0),
	mBigTheta(new Matrix(R.size()[0]*Q.size()[0]+Q.size()[0]*(Q.size()[0]+1)/2,R.size()[0]*Q.size()[0]+Q.size()[0]*(Q.size()[0]+1)/2)),
	mBigr(new Matrix(Q.size()[0]*Q.size()[0],R.size()[0]*Q.size()[0]+Q.size()[0]*(Q.size()[0]+1)/2)),
	mADPalg(nullptr)
{

	std::shared_ptr<AlgorithmADP> tem (new T());
	mADPalg = (tem->Creat(mQ,mR,mP,mK0,stepf));
	long double eps = 1e-10;
	mThetaInv=mThetaInv+1/eps; 
}

	template <typename T>
	ControllerADP<T>::ControllerADP(const SymmetricMatrix& Q, const SymmetricMatrix& R, const double delta, const SymmetricMatrix& P, const Matrix& K, Step* stepf)
:mK0(K), mKadp(K), mQ(Q), mR(R), mP(P), mdelta(delta), itx(mxx.begin()), itxx(mIxx.begin()), itxu(mIxu.begin()),
	mThetaInv(R.size()[0]*Q.size()[0]+Q.size()[0]*(Q.size()[0]+1)/2, R.size()[0]*Q.size()[0]+Q.size()[0]*(Q.size()[0]+1)/2), 
	mBigV(R.size()[0]*Q.size()[0]+Q.size()[0]*(Q.size()[0]+1)/2,0),
	mBigTheta(new Matrix(R.size()[0]*Q.size()[0]+Q.size()[0]*(Q.size()[0]+1)/2,R.size()[0]*Q.size()[0]+Q.size()[0]*(Q.size()[0]+1)/2)),
	mBigr(new Matrix(Q.size()[0]*Q.size()[0],R.size()[0]*Q.size()[0]+Q.size()[0]*(Q.size()[0]+1)/2)),
	mADPalg(nullptr)
{
	std::shared_ptr<AlgorithmADP> tem (new T());
	mADPalg = (tem->Creat(mQ,mR,mP,mK0,stepf));
	long double eps = 1e-10;
	mThetaInv=mThetaInv+1/eps; 
}

	template <>
const std::vector<double> ControllerADP<AlgorithmPI>::input(const std::vector<double>& x, const double dt,  const double t, noise noif)
{	

	if (dt<1e-20) 
		return noif(mR.size()[0],t);



	if(mxx.size()==0) 
	{
		mxx.push_back(vecs(prod(x,x)));
		++itx;
	}else {
		mxx.push_back(vecs(prod(x,x)));
	}

	if(mIxx.size()==0) 
	{
		mIxx.push_back(-dt*vec(kProd(x,x)));
		//std::cout << "test post" << std::endl;
		++itxx;
	}
	else {
		//std::vector<double> tem(mIxx.back()+dt*vecs(prod(x,x)));
		//mIxx.push_back(tem);
		mIxx.push_back(mIxx.back()-dt*vec(kProd(x,x)));
	}


	//Matrix matOut(mK0 * x);
	//std::vector<double> u(mR.size()[0],sin(t));
	//noise noif=&sinusoidal;
	std::vector<double> u = noif(mR.size()[0],t);
	//std::vector<double> u(vec(mK0 * x));

	//mu.push_back(u);

	if(mIxu.size()==0) 
	{
		mIxu.push_back(-2*dt*vec(kProd(x,mR*u)));
		++itxu;
	} else{
		//std::vector<double> tem2=vec(mIxu.back()+2*dt*kProd(x,mR*u));
		//mIxu.push_back(tem2);
		mIxu.push_back(mIxu.back()-2*dt*vec(kProd(x,mR*u)));
	}

	//std::cout << mxx.size() << ',' <<t << std::endl;


	if (t > mdelta && mxx.size()>1 && mIxx.size()>1 && mIxu.size()>1){
		//std::cout << mdelta << std::endl;
		mdelta +=0.1;
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
		std::vector<double> vec0 = mxx.back() - *itx;
		std::vector<double> vec1 = mIxx.back() - *itxx;
		std::vector<double> vec2 = mIxu.back() - *itxu;

		std::vector<Matrix> optResult = mADPalg->online(vec0, vec1, vec2, mBigr, mThetaInv, mBigV);
		//double err= optResult[0].F();
		mP = optResult[1];
		mKadp= optResult[2];
		disp(mP);
		//mBigV = vec(mThetaInv * *mBigr * vec(mP));

		//std::cout << "the error is " << err << std::endl;
		//mP.disp();


		itx = mxx.end();
		itxx = mIxx.end();
		itxu = mIxu.end();
		--itx;
		--itxx;
		--itxu;

		//}

	}

	return u;

}
	template <typename T>
const std::vector<double> ControllerADP<T>::input(const std::vector<double>& x, const double dt,  const double t, noise noif)
{	
	if(mxx.size()==0) 
	{
		mxx.push_back(vec(kProd(x,x)));
		++itx;
	}else {
		mxx.push_back(vec(kProd(x,x)));
	}

	if(mIxx.size()==0) 
	{
		mIxx.push_back(dt*vecs(prod(x,x)));
		//std::cout << "test post" << std::endl;
		++itxx;
	}
	else {
		//std::vector<double> tem(mIxx.back()+dt*vecs(prod(x,x)));
		//mIxx.push_back(tem);
		mIxx.push_back(mIxx.back()+dt*vecs(prod(x,x)));
	}


	//Matrix matOut(mK0 * x);
	//std::vector<double> u(mR.size()[0],sin(t));
	//noise noif=&sinusoidal;
	std::vector<double> u = noif(mR.size()[0],t);
	//std::vector<double> u(vec(mK0 * x));

	//mu.push_back(u);

	if(mIxu.size()==0) 
	{
		mIxu.push_back(2*dt*vec(kProd(x,mR*u)));
		++itxu;
	} else{
		//std::vector<double> tem2=vec(mIxu.back()+2*dt*kProd(x,mR*u));
		//mIxu.push_back(tem2);
		mIxu.push_back(mIxu.back()+2*dt*vec(kProd(x,mR*u)));
	}

	//std::cout << mxx.size() << ',' <<t << std::endl;


	if (t > mdelta && mxx.size()>1 && mIxx.size()>1 && mIxu.size()>1){
		//std::cout << mdelta << std::endl;
		mdelta +=0.1;
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




		// this is the new one
		std::vector<double> vec0 = mxx.back() - *itx;
		std::vector<double> vec1 = mIxx.back() - *itxx;
		std::vector<double> vec2 = mIxu.back() - *itxu;

		std::vector<Matrix> optResult = mADPalg->online(vec0, vec1, vec2, mBigr, mThetaInv, mBigV);
		double err= optResult[0].F();
		mP = optResult[1];
		mKadp= optResult[2];
		mBigV = vec(mThetaInv * *mBigr * vec(mP));

		std::cout << "the error is " << err << std::endl;
		mP.disp();


		//if((mBigV-mBigVold)*(mBigV-mBigVold)<1e-10)
		if(err<1e-10)
		{
			long double eps = 1e-10;
			mThetaInv=mThetaInv*0+1/eps; 
			mBigV = mBigV * 0;
			mxx.clear();
			mIxx.clear();
			mIxu.clear();
			itx = mxx.begin();
			itxx = mIxx.begin();
			itxu = mIxu.begin();
			unsigned int n = mQ.size()[0];
			unsigned int m = mR.size()[0];
			mBigr.reset(new Matrix(n*n,m*n+n*(n+1)/2));
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

	return u;

}

template <typename T>
void ControllerADP<T>::dispAll(){

	std::cout <<"what is in mxx" << std::endl;
	for (auto itx = mxx.begin();itx!=mxx.end();++itx)
	{
		disp(*itx);
	}
	std::cout <<"what is in mIxx" << std::endl;
	for (auto itxx = mIxx.begin();itxx!=mIxx.end();++itxx)
	{
		disp(*itxx);
	}
	std::cout <<"what is in mIxu" << std::endl;
	for (auto itxu = mIxu.begin();itxu!=mIxu.end();++itxu)
	{
		disp(*itxu);
	}

}


	//template <typename T>
//void ControllerADP<T>::LS(const std::vector<double>& phi, const double d)
//{
	//mThetaInv = mThetaInv - 1 / (1 + double(t(phi)*mThetaInv*phi)) * mThetaInv * phi * t(mThetaInv * phi);
	///[>mBigr= *mBigr + phi * d;// half online case;
	//mBigV = mBigV +(d - mBigV * phi)* (vec(mThetaInv * phi));
	////std::cout<< double(mBigV*phi)-d << std::endl;
//}
