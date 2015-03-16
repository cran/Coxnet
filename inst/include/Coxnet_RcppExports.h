// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#ifndef __Coxnet_RcppExports_h__
#define __Coxnet_RcppExports_h__

#include <RcppEigen.h>
#include <Rcpp.h>

namespace Coxnet {

    using namespace Rcpp;

    namespace {
        void validateSignature(const char* sig) {
            Rcpp::Function require = Rcpp::Environment::base_env()["require"];
            require("Coxnet", Rcpp::Named("quietly") = true);
            typedef int(*Ptr_validate)(const char*);
            static Ptr_validate p_validate = (Ptr_validate)
                R_GetCCallable("Coxnet", "Coxnet_RcppExport_validate");
            if (!p_validate(sig)) {
                throw Rcpp::function_not_exported(
                    "C++ function with signature '" + std::string(sig) + "' not found in Coxnet");
            }
        }
    }

    inline List scaleC(Eigen::MatrixXd X) {
        typedef SEXP(*Ptr_scaleC)(SEXP);
        static Ptr_scaleC p_scaleC = NULL;
        if (p_scaleC == NULL) {
            validateSignature("List(*scaleC)(Eigen::MatrixXd)");
            p_scaleC = (Ptr_scaleC)R_GetCCallable("Coxnet", "Coxnet_scaleC");
        }
        RObject __result;
        {
            RNGScope __rngScope;
            __result = p_scaleC(Rcpp::wrap(X));
        }
        if (__result.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (__result.inherits("try-error"))
            throw Rcpp::exception(as<std::string>(__result).c_str());
        return Rcpp::as<List >(__result);
    }

    inline List OmegaC(Eigen::MatrixXd& Omega, Eigen::VectorXi& sgn) {
        typedef SEXP(*Ptr_OmegaC)(SEXP,SEXP);
        static Ptr_OmegaC p_OmegaC = NULL;
        if (p_OmegaC == NULL) {
            validateSignature("List(*OmegaC)(Eigen::MatrixXd&,Eigen::VectorXi&)");
            p_OmegaC = (Ptr_OmegaC)R_GetCCallable("Coxnet", "Coxnet_OmegaC");
        }
        RObject __result;
        {
            RNGScope __rngScope;
            __result = p_OmegaC(Rcpp::wrap(Omega), Rcpp::wrap(sgn));
        }
        if (__result.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (__result.inherits("try-error"))
            throw Rcpp::exception(as<std::string>(__result).c_str());
        return Rcpp::as<List >(__result);
    }

    inline List OmegaSC(Eigen::SparseMatrix<double>& OmegaS, Eigen::VectorXi& sgn) {
        typedef SEXP(*Ptr_OmegaSC)(SEXP,SEXP);
        static Ptr_OmegaSC p_OmegaSC = NULL;
        if (p_OmegaSC == NULL) {
            validateSignature("List(*OmegaSC)(Eigen::SparseMatrix<double>&,Eigen::VectorXi&)");
            p_OmegaSC = (Ptr_OmegaSC)R_GetCCallable("Coxnet", "Coxnet_OmegaSC");
        }
        RObject __result;
        {
            RNGScope __rngScope;
            __result = p_OmegaSC(Rcpp::wrap(OmegaS), Rcpp::wrap(sgn));
        }
        if (__result.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (__result.inherits("try-error"))
            throw Rcpp::exception(as<std::string>(__result).c_str());
        return Rcpp::as<List >(__result);
    }

    inline double max_lambdaC(Eigen::MatrixXd X, Eigen::VectorXd tevent, int N, Eigen::VectorXi nevent, Eigen::VectorXi nevent1, Eigen::VectorXi loc1, int n, double alpha, Eigen::VectorXd wbeta, int N0) {
        typedef SEXP(*Ptr_max_lambdaC)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_max_lambdaC p_max_lambdaC = NULL;
        if (p_max_lambdaC == NULL) {
            validateSignature("double(*max_lambdaC)(Eigen::MatrixXd,Eigen::VectorXd,int,Eigen::VectorXi,Eigen::VectorXi,Eigen::VectorXi,int,double,Eigen::VectorXd,int)");
            p_max_lambdaC = (Ptr_max_lambdaC)R_GetCCallable("Coxnet", "Coxnet_max_lambdaC");
        }
        RObject __result;
        {
            RNGScope __rngScope;
            __result = p_max_lambdaC(Rcpp::wrap(X), Rcpp::wrap(tevent), Rcpp::wrap(N), Rcpp::wrap(nevent), Rcpp::wrap(nevent1), Rcpp::wrap(loc1), Rcpp::wrap(n), Rcpp::wrap(alpha), Rcpp::wrap(wbeta), Rcpp::wrap(N0));
        }
        if (__result.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (__result.inherits("try-error"))
            throw Rcpp::exception(as<std::string>(__result).c_str());
        return Rcpp::as<double >(__result);
    }

    inline double pletaCm(Eigen::VectorXd& xb, Eigen::VectorXd& exb, Eigen::VectorXi& nevent, Eigen::VectorXi& nevent1, Eigen::VectorXi& loc1, int& n, int& ifast, int& itwo) {
        typedef SEXP(*Ptr_pletaCm)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_pletaCm p_pletaCm = NULL;
        if (p_pletaCm == NULL) {
            validateSignature("double(*pletaCm)(Eigen::VectorXd&,Eigen::VectorXd&,Eigen::VectorXi&,Eigen::VectorXi&,Eigen::VectorXi&,int&,int&,int&)");
            p_pletaCm = (Ptr_pletaCm)R_GetCCallable("Coxnet", "Coxnet_pletaCm");
        }
        RObject __result;
        {
            RNGScope __rngScope;
            __result = p_pletaCm(Rcpp::wrap(xb), Rcpp::wrap(exb), Rcpp::wrap(nevent), Rcpp::wrap(nevent1), Rcpp::wrap(loc1), Rcpp::wrap(n), Rcpp::wrap(ifast), Rcpp::wrap(itwo));
        }
        if (__result.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (__result.inherits("try-error"))
            throw Rcpp::exception(as<std::string>(__result).c_str());
        return Rcpp::as<double >(__result);
    }

    inline Eigen::VectorXd cvtrimC(Eigen::VectorXd beta, int nn, int nn2, Eigen::VectorXi loco, Eigen::MatrixXd XF, int NF, Eigen::VectorXi neventF, Eigen::VectorXi nevent1F, Eigen::VectorXi loc1F, int nF, Eigen::MatrixXd X, int N, Eigen::VectorXi nevent, Eigen::VectorXi nevent1, Eigen::VectorXi loc1, int n, int ifast, int itwo) {
        typedef SEXP(*Ptr_cvtrimC)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_cvtrimC p_cvtrimC = NULL;
        if (p_cvtrimC == NULL) {
            validateSignature("Eigen::VectorXd(*cvtrimC)(Eigen::VectorXd,int,int,Eigen::VectorXi,Eigen::MatrixXd,int,Eigen::VectorXi,Eigen::VectorXi,Eigen::VectorXi,int,Eigen::MatrixXd,int,Eigen::VectorXi,Eigen::VectorXi,Eigen::VectorXi,int,int,int)");
            p_cvtrimC = (Ptr_cvtrimC)R_GetCCallable("Coxnet", "Coxnet_cvtrimC");
        }
        RObject __result;
        {
            RNGScope __rngScope;
            __result = p_cvtrimC(Rcpp::wrap(beta), Rcpp::wrap(nn), Rcpp::wrap(nn2), Rcpp::wrap(loco), Rcpp::wrap(XF), Rcpp::wrap(NF), Rcpp::wrap(neventF), Rcpp::wrap(nevent1F), Rcpp::wrap(loc1F), Rcpp::wrap(nF), Rcpp::wrap(X), Rcpp::wrap(N), Rcpp::wrap(nevent), Rcpp::wrap(nevent1), Rcpp::wrap(loc1), Rcpp::wrap(n), Rcpp::wrap(ifast), Rcpp::wrap(itwo));
        }
        if (__result.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (__result.inherits("try-error"))
            throw Rcpp::exception(as<std::string>(__result).c_str());
        return Rcpp::as<Eigen::VectorXd >(__result);
    }

    inline List coxenetC(Eigen::MatrixXd X, Eigen::VectorXd tevent, double alpha, Eigen::VectorXd lambda, int nlambda, Eigen::VectorXd wbeta, int N, Eigen::VectorXi nevent, Eigen::VectorXi nevent1, Eigen::VectorXi loc1, int n, int p, int N0, double thresh, int maxit, int ifast) {
        typedef SEXP(*Ptr_coxenetC)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_coxenetC p_coxenetC = NULL;
        if (p_coxenetC == NULL) {
            validateSignature("List(*coxenetC)(Eigen::MatrixXd,Eigen::VectorXd,double,Eigen::VectorXd,int,Eigen::VectorXd,int,Eigen::VectorXi,Eigen::VectorXi,Eigen::VectorXi,int,int,int,double,int,int)");
            p_coxenetC = (Ptr_coxenetC)R_GetCCallable("Coxnet", "Coxnet_coxenetC");
        }
        RObject __result;
        {
            RNGScope __rngScope;
            __result = p_coxenetC(Rcpp::wrap(X), Rcpp::wrap(tevent), Rcpp::wrap(alpha), Rcpp::wrap(lambda), Rcpp::wrap(nlambda), Rcpp::wrap(wbeta), Rcpp::wrap(N), Rcpp::wrap(nevent), Rcpp::wrap(nevent1), Rcpp::wrap(loc1), Rcpp::wrap(n), Rcpp::wrap(p), Rcpp::wrap(N0), Rcpp::wrap(thresh), Rcpp::wrap(maxit), Rcpp::wrap(ifast));
        }
        if (__result.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (__result.inherits("try-error"))
            throw Rcpp::exception(as<std::string>(__result).c_str());
        return Rcpp::as<List >(__result);
    }

    inline List cvcoxenetC(Eigen::MatrixXd X, Eigen::VectorXd tevent, double alpha, Eigen::VectorXd lambda, int nlambda, Eigen::VectorXd wbeta, int N, Eigen::VectorXi nevent, Eigen::VectorXi nevent1, Eigen::VectorXi loc1, int n, int p, int N0, double thresh, int maxit, int ifast, Eigen::MatrixXd XF, int NF, Eigen::VectorXi neventF, Eigen::VectorXi nevent1F, Eigen::VectorXi loc1F, int nF) {
        typedef SEXP(*Ptr_cvcoxenetC)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_cvcoxenetC p_cvcoxenetC = NULL;
        if (p_cvcoxenetC == NULL) {
            validateSignature("List(*cvcoxenetC)(Eigen::MatrixXd,Eigen::VectorXd,double,Eigen::VectorXd,int,Eigen::VectorXd,int,Eigen::VectorXi,Eigen::VectorXi,Eigen::VectorXi,int,int,int,double,int,int,Eigen::MatrixXd,int,Eigen::VectorXi,Eigen::VectorXi,Eigen::VectorXi,int)");
            p_cvcoxenetC = (Ptr_cvcoxenetC)R_GetCCallable("Coxnet", "Coxnet_cvcoxenetC");
        }
        RObject __result;
        {
            RNGScope __rngScope;
            __result = p_cvcoxenetC(Rcpp::wrap(X), Rcpp::wrap(tevent), Rcpp::wrap(alpha), Rcpp::wrap(lambda), Rcpp::wrap(nlambda), Rcpp::wrap(wbeta), Rcpp::wrap(N), Rcpp::wrap(nevent), Rcpp::wrap(nevent1), Rcpp::wrap(loc1), Rcpp::wrap(n), Rcpp::wrap(p), Rcpp::wrap(N0), Rcpp::wrap(thresh), Rcpp::wrap(maxit), Rcpp::wrap(ifast), Rcpp::wrap(XF), Rcpp::wrap(NF), Rcpp::wrap(neventF), Rcpp::wrap(nevent1F), Rcpp::wrap(loc1F), Rcpp::wrap(nF));
        }
        if (__result.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (__result.inherits("try-error"))
            throw Rcpp::exception(as<std::string>(__result).c_str());
        return Rcpp::as<List >(__result);
    }

    inline List coxnetC(Eigen::MatrixXd& X, Eigen::VectorXd tevent, double alpha, Eigen::VectorXd lambda, int nlambda, Eigen::VectorXd wbeta, Eigen::SparseMatrix<double>& Omega, Eigen::MatrixXd loc, Eigen::VectorXi nadj, int N, Eigen::VectorXi nevent, Eigen::VectorXi nevent1, Eigen::VectorXi loc1, int n, int p, int N0, double thresh, int maxit, int ifast) {
        typedef SEXP(*Ptr_coxnetC)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_coxnetC p_coxnetC = NULL;
        if (p_coxnetC == NULL) {
            validateSignature("List(*coxnetC)(Eigen::MatrixXd&,Eigen::VectorXd,double,Eigen::VectorXd,int,Eigen::VectorXd,Eigen::SparseMatrix<double>&,Eigen::MatrixXd,Eigen::VectorXi,int,Eigen::VectorXi,Eigen::VectorXi,Eigen::VectorXi,int,int,int,double,int,int)");
            p_coxnetC = (Ptr_coxnetC)R_GetCCallable("Coxnet", "Coxnet_coxnetC");
        }
        RObject __result;
        {
            RNGScope __rngScope;
            __result = p_coxnetC(Rcpp::wrap(X), Rcpp::wrap(tevent), Rcpp::wrap(alpha), Rcpp::wrap(lambda), Rcpp::wrap(nlambda), Rcpp::wrap(wbeta), Rcpp::wrap(Omega), Rcpp::wrap(loc), Rcpp::wrap(nadj), Rcpp::wrap(N), Rcpp::wrap(nevent), Rcpp::wrap(nevent1), Rcpp::wrap(loc1), Rcpp::wrap(n), Rcpp::wrap(p), Rcpp::wrap(N0), Rcpp::wrap(thresh), Rcpp::wrap(maxit), Rcpp::wrap(ifast));
        }
        if (__result.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (__result.inherits("try-error"))
            throw Rcpp::exception(as<std::string>(__result).c_str());
        return Rcpp::as<List >(__result);
    }

    inline List cvcoxnetC(Eigen::MatrixXd& X, Eigen::VectorXd tevent, double alpha, Eigen::VectorXd lambda, int nlambda, Eigen::VectorXd wbeta, Eigen::SparseMatrix<double>& Omega, Eigen::MatrixXd loc, Eigen::VectorXi nadj, int N, Eigen::VectorXi nevent, Eigen::VectorXi nevent1, Eigen::VectorXi loc1, int n, int p, int N0, double thresh, int maxit, int ifast, Eigen::MatrixXd XF, int NF, Eigen::VectorXi neventF, Eigen::VectorXi nevent1F, Eigen::VectorXi loc1F, int nF) {
        typedef SEXP(*Ptr_cvcoxnetC)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_cvcoxnetC p_cvcoxnetC = NULL;
        if (p_cvcoxnetC == NULL) {
            validateSignature("List(*cvcoxnetC)(Eigen::MatrixXd&,Eigen::VectorXd,double,Eigen::VectorXd,int,Eigen::VectorXd,Eigen::SparseMatrix<double>&,Eigen::MatrixXd,Eigen::VectorXi,int,Eigen::VectorXi,Eigen::VectorXi,Eigen::VectorXi,int,int,int,double,int,int,Eigen::MatrixXd,int,Eigen::VectorXi,Eigen::VectorXi,Eigen::VectorXi,int)");
            p_cvcoxnetC = (Ptr_cvcoxnetC)R_GetCCallable("Coxnet", "Coxnet_cvcoxnetC");
        }
        RObject __result;
        {
            RNGScope __rngScope;
            __result = p_cvcoxnetC(Rcpp::wrap(X), Rcpp::wrap(tevent), Rcpp::wrap(alpha), Rcpp::wrap(lambda), Rcpp::wrap(nlambda), Rcpp::wrap(wbeta), Rcpp::wrap(Omega), Rcpp::wrap(loc), Rcpp::wrap(nadj), Rcpp::wrap(N), Rcpp::wrap(nevent), Rcpp::wrap(nevent1), Rcpp::wrap(loc1), Rcpp::wrap(n), Rcpp::wrap(p), Rcpp::wrap(N0), Rcpp::wrap(thresh), Rcpp::wrap(maxit), Rcpp::wrap(ifast), Rcpp::wrap(XF), Rcpp::wrap(NF), Rcpp::wrap(neventF), Rcpp::wrap(nevent1F), Rcpp::wrap(loc1F), Rcpp::wrap(nF));
        }
        if (__result.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (__result.inherits("try-error"))
            throw Rcpp::exception(as<std::string>(__result).c_str());
        return Rcpp::as<List >(__result);
    }

}

#endif // __Coxnet_RcppExports_h__