import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from cvxopt import solvers, matrix
solvers.options['show_progress'] = False


class SyntheticControl:
    """
    The synthetic control method is a statistical method used to evaluate the effect
    of an intervention in comparative case studies. It involves the construction of 
    a weighted combination of groups used as controls, to which the treatment group is 
    compared.


    Parameters
    ----------
    outer_loop_method : string, default = 'Nelder-Mead'
        Optimisation procedure used in outer scipy.minimize loop. Options are [‘Nelder-Mead’, 
        ‘Powell’, ‘CG’, ‘BFGS’, ‘Newton-CG’, ‘L-BFGS-B’, ‘TNC’, ‘COBYLA’, ‘SLSQP’, 
        ‘trust-constr’, ‘dogleg’, ‘trust-ncg’, ‘trust-exact’, ‘trust-krylov’, custom] 
        
    maxiter : int, default = 100000
        Maximum number of iterations in optimisation.
        
    maxfev : int, default = 100000
        Maximum number of function evaluations in optimisation.
        
    xatol : float, default = 1e-3
        X tolerance for optimisation loop.
        
    fatol : float, default = 1e-12
        Function tolerance for optimisation loop.
    
    verbose_outer : bool, default = True
        Boolean to print or suppress epochs in outer loop.
    
    verbose_inner : bool, default = True
        Boolean to print or suppress epochs in inner loop.
        
    weight_nonnegative_restr : bool, default = True
        Whether to include nonnegative restriction on weights.
        
    weight_sum_restr : bool, default = True
        Whether to include restriction of summing to 1 of weights.
    
        
    Attributes
    ----------
    n_donorcountries : int
        Number of donor countries seen during method `fit`.
        
    success : bool
        Whether or not the (outer) optimizer exited successfully.
    
    message : string
        Description of the cause of the termination of (outer) optimization.
    
    weights : ndarray of shape n_donorcountries
        Fitted coefficient of the features in the decision function. When fit_intercept is
        set to true, coefficients includes bias term, hence first value is intercept.
        
    
    Methods
    ----------
    loss_function(v2, y, X)
        Loss function for inner loop.
    
    fit(X,y)
        Fit Synthetic Control Method and estimate weights.
        
    predict_doppelganger(X)
        Predict target variable for doppelganger using the fitted donor country weights 
        of SCM estimator.

        
    Examples
    --------
    >>> import pandas as pd
    >>> import numpy as np
    >>> from scipy.optimize import minimize
    >>> from cvxopt import solvers, matrix
    >>> solvers.options['show_progress'] = False

    >>> X = np.array(pd.read_excel('X0.xlsx', header=None))
    >>> y = np.array(pd.read_excel('X1.xlsx', header=None))

    >>> SCfit = SyntheticControl().fit(X,y)
    >>> np.round(SCfit.weights,2)
    array([[0.  , 0.  , 0.  , 0.13, 0.  , 0.  , 
            0.  , 0.  , 0.  , 0.  , 0.08, 0.  , 
            0.1 , 0.08, 0.  , 0.  , 0.07, 0.  , 
            0.19, 0.07, 0.  , 0.  , 0.  , 0.29]])
    """
    
    def __init__(
        self,
        outer_loop_method = 'Nelder-Mead',
        maxiter = 100000,
        maxfev = 100000,
        xatol = 1e-3,
        fatol = 1e-12,
        verbose_outer = True,
        verbose_inner = True,
        weight_nonnegative_restr = True,
        weight_sum_restr = True,
        count = 0,
    ):
        self.outer_loop_method = outer_loop_method
        self.maxiter = maxiter
        self.maxfev = maxfev
        self.xatol = xatol
        self.fatol = fatol
        self.verbose_outer = verbose_outer
        self.verbose_inner = verbose_inner
        self.weight_nonnegative_restr = weight_nonnegative_restr
        self.weight_sum_restr = weight_sum_restr
        self.count = count
        
    def loss_function(self, v2, y, X):
        """
        Define loss function for inner loop.

        Parameters
        ----------
        v2 : array of shape (n_samples - 1)
            Initial values for V matrix.
            
        y : array of shape (n_samples)
            Data of treatment country, stacked in vector.
            
        X : array of shape (n_samples, n_donorcountries)
            Data of donor countries, stacked in matrix.
            

        Returns
        -------
        SSR : float
            Sum of squared residuals after optimization inner loop.
        """
        
        # Initalize matrices
        V = np.diag(np.insert(np.exp(v2), 0, 1))
        H = np.transpose(X) @ V @ X
        l = np.shape(X)[1]
        
        # Create correct CVxOPT matrices for quad. optimization loop with constraints
        Q = matrix((H + np.transpose(H))/2)
        r = matrix((-np.transpose(y) @ V @ X).T)
        G = matrix(-np.eye(l))
        h = matrix(np.zeros(l))
        A = matrix(np.ones(l)).T 
        b = matrix(1.0)

        # Optimization
        if self.weight_nonnegative_restr == False:
            G, h  = None, None
        if self.weight_sum_restr == False:
            A, b = None, None
            
        w = solvers.qp(Q, r, G, h, A, b)['x']
        SSR = np.sum((y - X@w)**2)

        # Print progress
        self.count += 1
        if (self.verbose_inner == True) & (self.count%1000==0):
            print(f'Function evaluation {self.count}/{self.maxfev} \tloss = {SSR}')
            
        return SSR
        
    def fit(self, X, y):
        """
        Fit Synthetic Control Method.

        Parameters
        ----------
        X : array of shape (n_samples, n_donorcountries)
            Data of donor countries, stacked in matrix.
        y : array of shape (n_samples)
            Data of treatment country, stacked in vector.

        Returns
        -------
        self : object
            Fitted Estimator.
        """
        
        # Find number of donor countries
        self.n_donorcountries = np.shape(X)[1]
        
        # Set starting values for outer loop
        s = np.hstack((y, X)).std(axis=1)
        s1, s2 = s[0], s[1:]
        v20 = np.log((s1/s2)**2)
        minimize_options = {'maxiter': self.maxiter, 
                            'maxfev': self.maxfev,
                            'disp': self.verbose_outer,
                            'xatol': self.xatol, 
                            'fatol': self.fatol}
        
        # Outer loop
        result = minimize(self.loss_function, v20, args=(y, X), method=self.outer_loop_method, 
              options = minimize_options)
        v2 = result.x
        
        # Store results outer loop
        self.success = result.success
        self.message = result.message
        
        # Restore weights from V weights
        self.V = np.diag(np.insert(np.exp(v2), 0, 1))
        H = np.transpose(X) @ self.V @ X
        l = np.shape(X)[1]
        
        # Transform to correct CVxOPT types
        Q = matrix((H + np.transpose(H))/2)
        r = matrix((-np.transpose(y) @ self.V @ X).T)
        G = matrix(-np.eye(l))
        h = matrix(np.zeros(l))
        A = matrix(np.ones(l)).T   
        b = matrix(1.0)
        
        # Optimization and resulting weights
        if self.weight_nonnegative_restr == False:
            G, h  = None, None
        if self.weight_sum_restr == False:
            A, b = None, None
        
        self.weights = np.array(solvers.qp(Q, r, G, h, A, b)['x']).T

        return self
        
       
    def predict_doppelganger(self, X):
        """
        Predict target variable for doppelganger using the fitted donor country 
        weights of SCM estimator.

        Parameters
        ----------
        X : array of shape (n_samples, n_donorcountries)
            Values for features of new data.

        Returns
        -------
        y_predicted : array of shape (n_samples)
            Returns predicted values for donor countries.
        """
        
        if hasattr(self, 'weights') == False:
                raise ValueError(
                     " This SyntheticControl instance is not fitted yet." \
                     " Call 'fit' method with appropriate X and y before using this predict function."
                 )
        
        y_predicted = self.weights @ X.T
        
        return y_predicted.T
