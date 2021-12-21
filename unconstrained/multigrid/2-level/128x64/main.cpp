#include <iostream>
#include <fstream>
#include <eigen/Eigen/Dense>
#include <vector>
#include <math.h>

// Main functions starts from here
// function to read txt file and adding values to vector elements
Eigen::VectorXf CSVtoVector(std::string file, int i, int n) {

    double data[i*n];

    Eigen::VectorXf vector = Eigen::VectorXf::Zero(i*n);

    std::ifstream input(file);

    for (int j = 0; j < i*n; ++j) {
        input >> data[j];
        vector[j] = data[j];

    }
    
    // std::cout << vector << std::endl;

    return vector;
}

// function to read txt file and adding values to vector elements
Eigen::VectorXf write_w(std::string file, int i, int n) {

    double data[2*(n+1)*(i-1)];

    Eigen::VectorXf vector = Eigen::VectorXf::Zero(2*(n+1)*(i-1));

    std::ifstream input(file);

    for (int j = 0; j < 2*(n+1)*(i-1) - 1; ++j) {
        input >> data[j];
        vector[j] = data[j];

    }
    
    // std::cout << vector << std::endl;

    return vector;
}

// Saving the vector to the file
void WriteFile(Eigen::VectorXf vector, std::string filename) {
    
    std::ofstream myfile (filename.append(".txt"), std::ios_base::app);

    if (myfile.is_open()){

            myfile << vector << "\n";
    }
    myfile.close();
}

Eigen::VectorXf ProjectionFormula(Eigen::VectorXf p, int i, int n, float u_a, float u_b){

    for (int j = 0; j < (n+1)*(i-1); ++j) {
        
        if(p[j] > u_b){
            p[j] = u_b;
        }

        if(p[j] < u_a){
            p[j] = u_a;
        }
    }

    return p;
}

// Building Matrix A
Eigen::MatrixXf BuildA(int n, int i, float delt, float delx, float mu, Eigen::VectorXf u){

    Eigen::MatrixXf A(2*n*i + 2*i, 2*n*i + 2*i);

    // Eigen::VectorXf u = Eigen::VectorXf::Random((i-1)*(n+1));
    Eigen::MatrixXf u_matrix = u.asDiagonal();
    
    //Constructing matrix D
    Eigen::MatrixXf D =  Eigen::MatrixXf::Identity(i, i);

    for(int j = 0; j < i; ++j){

        D(j, j) = -2.0;

    }

    for(int j=0; j < i-1; ++j){

        D(j+1, j) = 1.0;

    }

    for(int j=0; j < i-1; ++j){

        D(j, j+1) = 1.0;

    }

    D = 1.0/(delx*delx) * D;

    // std::cout << "\n Matrix D: \n" << D << std::endl;
           
    // constructing submatrix A0 for initial time step
    Eigen::MatrixXf A0 = Eigen::MatrixXf::Zero(2*i, 2*i);

    Eigen::MatrixXf a01 = Eigen::MatrixXf::Identity(i, i);
    Eigen::MatrixXf a02 = Eigen::MatrixXf::Zero(i, i);
    Eigen::MatrixXf a11 = delt * Eigen::MatrixXf::Identity(i, i);
    Eigen::MatrixXf a12 = (delt * D) - (Eigen::MatrixXf::Identity(i, i)) + delt * u_matrix.block(0, 0, i, i);

    A0.block(0, 0, i, i) = a01;
    A0.block(0, i, i, i) = a02;
    A0.block(i, 0, i, i) = a11;
    A0.block(i, i, i, i) = a12; 

    for(int j=1; j<n; ++j){

        // constructing submatrix Am for all other time steps
        Eigen::MatrixXf Am = Eigen::MatrixXf::Zero(2*i, 2*i);

        Eigen::MatrixXf Am01 = (delt * D) - (Eigen::MatrixXf::Identity(i, i)) + delt * u_matrix.block(j*i, j*i, i, i);
        Eigen::MatrixXf Am02 = Eigen::MatrixXf::Zero(i, i);
        Eigen::MatrixXf Am11 = delt * Eigen::MatrixXf::Identity(i, i);
        Eigen::MatrixXf Am12 = (delt * D) - (Eigen::MatrixXf::Identity(i, i)) + delt * u_matrix.block(j*i, j*i, i, i);

        Am.block(0, 0, i, i) = Am01;
        Am.block(0, i, i, i) = Am02;
        Am.block(i, 0, i, i) = Am11;
        Am.block(i, i, i, i) = Am12;

        A.block(2*i*j, 2*i*j, 2*i, 2*i) = Am;
        
    }

    // constructing submatrix Ant for final time step
    Eigen::MatrixXf Ant = Eigen::MatrixXf::Zero(2*i, 2*i);

    Eigen::MatrixXf Ant01 = (delt * D) - (Eigen::MatrixXf::Identity(i, i)) + delt * u_matrix.block(n*i, n*i, i, i);
    Eigen::MatrixXf Ant02 = Eigen::MatrixXf::Zero(i, i);
    Eigen::MatrixXf Ant11 = Eigen::MatrixXf::Zero(i, i);
    Eigen::MatrixXf Ant12 = Eigen::MatrixXf::Identity(i, i);

    Ant.block(0, 0, i, i) = Ant01;
    Ant.block(0, i, i, i) = Ant02;
    Ant.block(i, 0, i, i) = Ant11;
    Ant.block(i, i, i, i) = Ant12;

    // constructing submatrix B
    Eigen::MatrixXf B = Eigen::MatrixXf::Zero(2*i, 2*i);
    Eigen::MatrixXf B01 = Eigen::MatrixXf::Identity(i, i);
    B.block(0, 0, i, i) = B01;

    // constructing submatrix C
    Eigen::MatrixXf C = Eigen::MatrixXf::Zero(2*i, 2*i);
    Eigen::MatrixXf C11 = Eigen::MatrixXf::Identity(i, i);
    C.block(i, i, i, i) = C11;

    // substituting all the submatrices to matrix A
    A.block(0, 0, 2*i, 2*i) = A0;
    A.block(2*n*i, 2*n*i, 2*i, 2*i) = Ant;

    for(int j=0; j < n; ++j){

        A.block(2*i*(1+j), 2*i*j, 2*i, 2*i) = B;
    }
        

    for(int j=0; j < n; ++j){

        A.block(2*i*j, 2*i*(1+j), 2*i, 2*i) = C;
    } 

    return A;
} 
// Building Matrix A ends

// Building Matrix M
Eigen::MatrixXf BuildM(int n,int i, float delt, float delx, float mu, Eigen::VectorXf u){

    Eigen::MatrixXf M(2*n*i + 2*i, 2*n*i + 2*i);

    //Eigen::VectorXf u = Eigen::VectorXf::Random((i-1)*(n+1));
    Eigen::MatrixXf u_matrix = u.asDiagonal();

    //Constructing matrix D
    Eigen::MatrixXf D =  Eigen::MatrixXf::Identity(i, i);

    for(int j = 0; j < i; ++j){

        D(j, j) = -2.0;

    }

    for(int j=0; j < i-1; ++j){

        D(j+1, j) = 1.0;

    }

    for(int j=0; j < i-1; ++j){

        D(j, j+1) = 1.0;

    }

    D = 1.0/(delx*delx) * D;
    
    // constructing submatrix A0 for initial time step
    Eigen::MatrixXf m0 = Eigen::MatrixXf::Zero(2*i, 2*i);

    Eigen::MatrixXf m01 = Eigen::MatrixXf::Identity(i, i);
    Eigen::MatrixXf m02 = Eigen::MatrixXf::Zero(i, i);
    Eigen::MatrixXf m11 = delt * Eigen::MatrixXf::Identity(i, i);
    Eigen::MatrixXf m12 = (delt * D) - (Eigen::MatrixXf::Identity(i, i)) + delt * u_matrix.block(0, 0, i, i);

    m0.block(0, 0, i, i) = m01;
    m0.block(0, i, i, i) = m02;
    m0.block(i, 0, i, i) = m11;
    m0.block(i, i, i, i) = m12; 

    for(int j=1; j<n; ++j){

        // constructing submatrix Am for all other time steps
        Eigen::MatrixXf m = Eigen::MatrixXf::Zero(2*i, 2*i);

        Eigen::MatrixXf m01 = (delt * D) - (Eigen::MatrixXf::Identity(i, i)) + delt * u_matrix.block(j*i, j*i, i, i);
        Eigen::MatrixXf m02 = Eigen::MatrixXf::Zero(i, i);
        Eigen::MatrixXf m11 = delt * Eigen::MatrixXf::Identity(i, i);
        Eigen::MatrixXf m12 = (delt * D) - (Eigen::MatrixXf::Identity(i, i)) + delt * u_matrix.block(j*i, j*i, i, i);

        m.block(0, 0, i, i) = m01;
        m.block(0, i, i, i) = m02;
        m.block(i, 0, i, i) = m11;
        m.block(i, i, i, i) = m12;

        M.block(2*i*j, 2*i*j, 2*i, 2*i) = m;
        
    }

    // constructing submatrix Ant for final time step
    Eigen::MatrixXf mnt = Eigen::MatrixXf::Zero(2*i, 2*i);

    Eigen::MatrixXf mnt01 = (delt * D) - (Eigen::MatrixXf::Identity(i, i)) + delt * u_matrix.block(n*i, n*i, i, i);
    Eigen::MatrixXf mnt02 = Eigen::MatrixXf::Zero(i, i);
    Eigen::MatrixXf mnt11 = Eigen::MatrixXf::Zero(i, i);
    Eigen::MatrixXf mnt12 = Eigen::MatrixXf::Identity(i, i);

    mnt.block(0, 0, i, i) = mnt01;
    mnt.block(0, i, i, i) = mnt02;
    mnt.block(i, 0, i, i) = mnt11;
    mnt.block(i, i, i, i) = mnt12;

    M.block(0, 0, 2*i, 2*i) = m0;
    M.block(2*n*i, 2*n*i, 2*i, 2*i) = mnt;

    return M;
} 
// Building Matrix M ends


// Building vector f
Eigen::VectorXf BoundaryConditions(Eigen::VectorXf f, Eigen::VectorXf Y, Eigen::VectorXf P, Eigen::VectorXf Z, Eigen::VectorXf f_i, int i, int n){
    
    float delt = 1.0/n;

    for(int j=0; j < n+1; ++j){

        f.block(2*i*j, 0, i, 1) =  delt * f_i.block(i*j, 0, i, 1);

        f.block(2*i*j + i, 0, i, 1) =   delt * Z.block(i*j, 0, i, 1);
        
    }
    
    f.block(0, 0, i, 1) = Y.block(0, 0, i, 1);
    f.block(2*i*n + i, 0, i, 1) = P.block(i*n, 0, i, 1);

    return f;
}
// Building vector f ends


// Performing LDU decomposition
std::vector<Eigen::MatrixXf> LUdecomposition(Eigen::MatrixXf M){

    std::vector<Eigen::MatrixXf> myVec;

    Eigen::MatrixXf U = Eigen::MatrixXf::Zero(M.rows(), M.cols());
    U.triangularView<Eigen::StrictlyUpper>() = M.triangularView<Eigen::StrictlyUpper>();
    
    Eigen::MatrixXf L = Eigen::MatrixXf::Zero(M.rows(), M.cols());
    L.triangularView<Eigen::StrictlyLower>() = M.triangularView<Eigen::StrictlyLower>();
    
    Eigen::MatrixXf D = Eigen::MatrixXf::Zero(M.rows(), M.cols());
    for(int i=0; i < M.rows(); ++i){

        D(i, i) = M(i, i);

    }
    
    myVec.push_back(L);
    myVec.push_back(D);
    myVec.push_back(U);

    return myVec;
}
// Performing LDU decomposition ends


// Building Jacobi Smoother
// Solving in sequential manner for each time step
Eigen::VectorXf SpaceTimeBlockJacobiSmoother(Eigen::MatrixXf A, Eigen::VectorXf w_old, Eigen::VectorXf f, int n, int i, int Iterations, float Tol){

    Eigen::VectorXf d;
    Eigen::VectorXf w_new;
    Eigen::VectorXf w_diff;

    float norm = 1.0;
    float omega = 0.01;
    
    std::vector<Eigen::MatrixXf> myVec= LUdecomposition(A);
    Eigen::MatrixXf L = myVec[0];
    Eigen::MatrixXf D = myVec[1];
    Eigen::MatrixXf U = myVec[2];

    // adding initial line to the file
    std::ofstream myfile ("Residual.txt", std::ios_base::app);

        if (myfile.is_open()){

            myfile << "Residual: " << "\n";
        }
        myfile.close();

    // adding initial line to the file
    std::ofstream myfile_2 ("Iterations.txt", std::ios_base::app);

        if (myfile_2.is_open()){

            myfile_2 << "Iterations: " << "\n";
        }
        myfile_2.close();    

    int j = 0;
    while( (norm > Tol) && (j < Iterations) ){
        
        // calculating the defect
        d = f - A * w_old;

        for(int k=0; k < n; ++k){
            
            d.block(2*i*k, 0, 2*i, 1) = D.block(2*i*k, 2*i*k, 2*i, 2*i).inverse() * d.block(2*i*k, 0, 2*i, 1);
              
        }

        // Update the vector
        w_new = w_old + omega * d;

        // calculating the norm
        w_diff = w_new - w_old;

        //std::cout << "\n w_diff: \n" << w_diff;

        norm = w_diff.norm();

        // Printing the residual norm
        if(j==0){
            std::cout << "\n Residual norm: " << "\n" <<  std::endl;
        }
        std::cout << norm << "\n" <<  std::endl;
        
        // writing residual norm to the file
        std::ofstream myfile ("Residual.txt", std::ios_base::app);

        if (myfile.is_open()){

            myfile << norm << "\n";
        }
        myfile.close();

        // writing number of iterations to the file
        std::ofstream myfile_2 ("Iterations.txt", std::ios_base::app);

        if (myfile_2.is_open()){

            myfile_2 << j << "\n";
        }
        myfile_2.close();

        w_old = w_new;

        ++j;
    }

    std::cout << "\n Solution Converged in: "  << j << " Iterations." << std::endl;

    return w_new;
}
// Building Jacobi Smoother ends


// Building Gauss Seidel Smoother
Eigen::VectorXf SpaceTimeBlockGaussSeidelSmoother(Eigen::MatrixXf A, Eigen::VectorXf w_old, Eigen::VectorXf f, int n, int i, int Iterations, float Tol){

    Eigen::VectorXf d;
    Eigen::VectorXf w_new;
    Eigen::VectorXf w_diff;

    float norm = 1;
    float omega = 0.01;

    std::vector<Eigen::MatrixXf> myVec= LUdecomposition(A);
    Eigen::MatrixXf L = myVec[0];
    Eigen::MatrixXf D = myVec[1];
    Eigen::MatrixXf U = myVec[2];

    // adding initial line to the file
    std::ofstream myfile ("Residual.txt", std::ios_base::app);

        if (myfile.is_open()){

            myfile << "Residual: " << "\n";
        }
        myfile.close();

    // adding initial line to the file
    std::ofstream myfile_2 ("Iterations.txt", std::ios_base::app);

        if (myfile_2.is_open()){

            myfile_2 << "Iterations: " << "\n";
        }
        myfile_2.close();    

    int j = 0;
    while( (j < Iterations) && (norm > Tol) ){
        
        // calculating the defect
        d = f - A * w_old;

        for(int k=0; k < n; ++k){
            
            d.block(2*i*k, 0, 2*i, 1) = ( L.block(2*i*k, 2*i*k, 2*i, 2*i) + D.block(2*i*k, 2*i*k, 2*i, 2*i)).inverse() * d.block(2*i*k, 0, 2*i, 1);
              
        }

        w_new = w_old + omega * d;
        
        // calculating the norm
        w_diff = w_new - w_old;

        norm = w_diff.norm();

        // Printing the norm
        if(j==0){
            std::cout << "\n Residual norm: " << "\n" <<  std::endl;
        }
        std::cout << norm << "\n" <<  std::endl;

        // Writing the residual norm to the file 
        std::ofstream myfile ("Residual.txt", std::ios_base::app);

        if (myfile.is_open()){

            myfile << norm << "\n";
        }
        myfile.close();

        // Writing number of iterations to the file 
        std::ofstream myfile_2 ("Iterations.txt", std::ios_base::app);

        if (myfile_2.is_open()){

            myfile_2 << j << "\n";
        }
        myfile_2.close();

        w_old = w_new;

        ++j;
    }
    
    std::cout << "\n Solution Converged in: "  << j << " Iterations " << std::endl;

    return w_new;

}
// Building Gauss Seidel Smoother ends

// Building SOR Smoother
Eigen::VectorXf SpaceTimeBlockSORSmoother(Eigen::MatrixXf A, Eigen::VectorXf w_old, Eigen::VectorXf f, float omega_2, int n, int i, int Iterations, float Tol){

    std::cout << "n :" << n << std::endl;
    std::cout << "i :" << i << std::endl;

    Eigen::VectorXf d;
    Eigen::VectorXf w_new;
    Eigen::VectorXf w_diff;

    float norm = 1;
    float omega_1 = 0.03;
    
    std::vector<Eigen::MatrixXf> myVec= LUdecomposition(A);
    Eigen::MatrixXf L = myVec[0];
    Eigen::MatrixXf D = myVec[1];
    Eigen::MatrixXf U = myVec[2];

    // adding initial line to the file
    std::ofstream myfile ("Residual.txt", std::ios_base::app);

        if (myfile.is_open()){

            myfile << "Residual: " << "\n";
        }
        myfile.close();

    // adding initial line to the file
    std::ofstream myfile_2 ("Iterations.txt", std::ios_base::app);

        if (myfile_2.is_open()){

            myfile_2 << "Iterations: " << "\n";
        }
        myfile_2.close();    

    int j = 0;
    while( (j < Iterations) && (norm > Tol) ){
        
        // calculating the defect
        d = f - A * w_old;

        for(int k=0; k < n; ++k){
            
            d.block(2*i*k, 0, 2*i, 1) = ( L.block(2*i*k, 2*i*k, 2*i, 2*i) + (1/omega_2) * D.block(2*i*k, 2*i*k, 2*i, 2*i)).inverse() * d.block(2*i*k, 0, 2*i, 1);
              
        }

        w_new = w_old + omega_1 * d;
        
        // calculating the norm
        w_diff = w_new - w_old;

        norm = w_diff.norm();

        // Printing the norm
        if(j==0){
            std::cout << "\n Residual norm: " << "\n" <<  std::endl;
        }
        std::cout << norm << "\n" <<  std::endl;

        // Writing the residual norm to the file 
        std::ofstream myfile ("Residual.txt", std::ios_base::app);

        if (myfile.is_open()){

            myfile << norm << "\n";
        }
        myfile.close();

        // Writing number of iterations to the file 
        std::ofstream myfile_2 ("Iterations.txt", std::ios_base::app);

        if (myfile_2.is_open()){

            myfile_2 << j << "\n";
        }
        myfile_2.close();

        w_old = w_new;

        ++j;
    }
    
    // std::cout << "\n Solution Converged in: "  << j << " Iterations " << std::endl;

    return w_new;

}
// Building SOR Smoother ends

// // Building SOR Smoother
Eigen::VectorXf SpaceTimeForwardBackwardBlockSORSmoother(Eigen::MatrixXf A, Eigen::VectorXf w_old, Eigen::VectorXf f, int n, int i, int NSM, float Tol){

    // std::cout << "n :" << n << std::endl;
    // std::cout << "i :" << i << std::endl;

    Eigen::VectorXf r;
    Eigen::VectorXf w_new;
    Eigen::VectorXf w_diff;

    float norm = 1;
    float omega_1 = 0.8;
    float omega_2 = 0.2;
    
    Eigen::VectorXf d = Eigen::VectorXf::Zero(2*i*(n+1));

    // adding initial line to the file
    std::ofstream myfile ("Residual.txt", std::ios_base::app);

        if (myfile.is_open()){

            myfile << "Residual: " << "\n";
        }
        myfile.close();

    // adding initial line to the file
    std::ofstream myfile_2 ("Iterations.txt", std::ios_base::app);

        if (myfile_2.is_open()){

            myfile_2 << "Iterations: " << "\n";
        }
        myfile_2.close();    

    Eigen::MatrixXf M_tilda = Eigen::MatrixXf::Zero(2*i, 2*i);
    Eigen::MatrixXf M_tilda_01 = Eigen::MatrixXf::Identity(i, i);
    M_tilda.block(0, 0, i, i) = M_tilda_01;

    Eigen::MatrixXf M_hat = Eigen::MatrixXf::Zero(2*i, 2*i);
    Eigen::MatrixXf M_hat_11 = Eigen::MatrixXf::Identity(i, i);
    M_hat.block(i, i, i, i) = M_hat_11;

    int j = 0;
    while( (j < NSM) && (norm > Tol) ){
        
        // calculating the defect
        r = f - A * w_old;

        Eigen::VectorXf x = Eigen::VectorXf::Zero(2*i*(n+1));
        Eigen::VectorXf x_old;

        for(int k=0; k < 1; k++){

            x_old = x;

            for(int l=0; l<=n; l++){
                
                Eigen::VectorXf _x, x_;                     // x_-1 = 0, x_(N+1) = 0

                if (l == 0){
                    _x = Eigen::VectorXf::Zero(2*i, 1);
                }
                else{
                    _x = x.block(2*i*(l-1), 0, 2*i, 1);
                }
                
                if (l == n){
                    x_ = Eigen::VectorXf::Zero(2*i, 1);
                }
                else{
                    x_ = x_old.block(2*i*(l+1), 0, 2*i, 1);
                }

                d.block(2*i*l, 0, 2*i, 1) = r.block(2*i*l, 0, 2*i, 1) - M_tilda * _x - M_hat * x_;
                x.block(2*i*l, 0, 2*i, 1) = (1-omega_1) * x_old.block(2*i*l, 0, 2*i, 1) + omega_1 * A.block(2*i*l, 2*i*l, 2*i, 2*i).inverse() * d.block(2*i*l, 0, 2*i, 1);
            }

            x_old = x;

            for(int l=n; l>=0; l--){
                
                Eigen::VectorXf _x, x_;                        // x_-1 = 0, x_(N+1) = 0

                if (l == 0){
                    _x = Eigen::VectorXf::Zero(2*i, 1);
                }
                else{
                    _x = x_old.block(2*i*(l-1), 0, 2*i, 1);
                }

                if (l == n){
                    x_ = Eigen::VectorXf::Zero(2*i, 1);
                }
                else{
                    x_ = x.block(2*i*(l+1), 0, 2*i, 1);
                }

                d.block(2*i*l, 0, 2*i, 1) = r.block(2*i*l, 0, 2*i, 1) - M_tilda * _x - M_hat * x_;
                x.block(2*i*l, 0, 2*i, 1) = (1-omega_1) * x_old.block(2*i*l, 0, 2*i, 1) + omega_1 * A.block(2*i*l, 2*i*l, 2*i, 2*i).inverse() * d.block(2*i*l, 0, 2*i, 1);
            }
              
        }
        
        w_new = w_old + omega_2 * x;
        
        // calculating the norm
        w_diff = w_new - w_old;

        norm = w_diff.norm();

        // Printing the norm
        if(j==0){
            std::cout << "\n Residual norm: " << "\n" <<  std::endl;
        }
        std::cout << norm << "\n" <<  std::endl;

        // Writing the residual norm to the file 
        std::ofstream myfile ("Residual.txt", std::ios_base::app);

        if (myfile.is_open()){

            myfile << norm << "\n";
        }
        myfile.close();

        // Writing number of iterations to the file 
        std::ofstream myfile_2 ("Iterations.txt", std::ios_base::app);

        if (myfile_2.is_open()){

            myfile_2 << j << "\n";
        }
        myfile_2.close();

        w_old = w_new;

        ++j;
    }

    //std::cout << "\n Solution Converged in: "  << j << " Iterations " << std::endl;

    return w_new;
}
// // Building SOR Smoother ends

// Restriction starts here
// Space Restriction operator starts here
Eigen::MatrixXf Restriction_Operator(int n, int i){

    Eigen::VectorXf d_space_restricted;
        
    Eigen::MatrixXf d_ = Eigen::MatrixXf::Zero((i/2 - 1), i-1);
    
    Eigen::MatrixXf d_m(1, 3);

    d_m << 1, 2, 1;
    
    for(int j=0; j < (i/2 - 1); ++j){

        d_.block(j, 2*j, 1, 3) = d_m;
    }

    Eigen::MatrixXf Operator = Eigen::MatrixXf::Zero((i/2 - 1) * 2, (i-1) * 2);

    Operator.block(0, 0, (i/2 - 1), i-1) = d_;
    Operator.block((i/2 - 1), i-1, (i/2 - 1), i-1) = d_;

    return Operator;
}
// Space Restriction operator ends here

// Space Restriction operator starts here
Eigen::VectorXf Space_Restriction(Eigen::VectorXf d, int i){

    Eigen::VectorXf d_space_restricted;
    
    int n = d.size()/((i-1)*2);

    Eigen::MatrixXf space_res_operator = Eigen::MatrixXf::Zero( (i/2 - 1) * 2 * n, d.size());
    
    Eigen::MatrixXf d_ = Eigen::MatrixXf::Zero((i/2 - 1), i-1);
    
    Eigen::MatrixXf d_m(1, 3);

    d_m << 1, 2, 1;
    
    for(int j=0; j < (i/2 - 1); ++j){

        d_.block(j, 2*j, 1, 3) = d_m;
    }

    Eigen::MatrixXf time_block = Eigen::MatrixXf::Zero((i/2 - 1) * 2, (i-1) * 2);

    time_block.block(0, 0, (i/2 - 1), i-1) = d_;
    time_block.block((i/2 - 1), i-1, (i/2 - 1), i-1) = d_;

    //std::cout << "\n restriction operator: \n" << time_block;

    for(int j=0; j < n; ++j){

        space_res_operator.block((i/2 - 1) * 2 * j, (i-1) * 2 * j, (i/2 - 1) * 2, (i-1) * 2) = time_block;
    }

    //std::cout << "\n space_res_operator: \n" << space_res_operator;

    d_space_restricted = (1.0/4) * space_res_operator * d;

    return d_space_restricted;
}
// Restriction operator ends here

// Time Restriction operator starts here
Eigen::VectorXf Time_Restriction(Eigen::VectorXf d, int n, int i){

    Eigen::VectorXf d_time_restricted;
      
    Eigen::MatrixXf time_res_operator = Eigen::MatrixXf::Zero(((n+1)/2 + 1) * 2 * (i-1), d.size());
    
    Eigen::MatrixXf a_0 = Eigen::MatrixXf::Identity((i-1) * 2, (i-1) * 2 * 2);

    a_0.block(0, 0, (i-1) * 2, (i-1) * 2) = 2 * Eigen::MatrixXf::Identity((i-1) * 2, (i-1) * 2);
    a_0.block(0, (i-1) * 2, (i-1) * 2, (i-1) * 2) = Eigen::MatrixXf::Identity((i-1) * 2, (i-1) * 2);

    Eigen::MatrixXf a_m = Eigen::MatrixXf::Zero((i-1) * 2, (i-1) * 2 * 3);

    a_m.block(0, 0, (i-1) * 2, (i-1) * 2) = Eigen::MatrixXf::Identity((i-1) * 2, (i-1) * 2);
    a_m.block(0, (i-1) * 2, (i-1) * 2, (i-1) * 2) = 2 * Eigen::MatrixXf::Identity((i-1) * 2, (i-1) * 2);
    a_m.block(0, (i-1) * 2 * 2, (i-1) * 2, (i-1) * 2) = Eigen::MatrixXf::Identity((i-1) * 2, (i-1) * 2);

    Eigen::MatrixXf a_nt = Eigen::MatrixXf::Identity((i-1) * 2, (i-1) * 2 * 2);

    a_nt.block(0, 0, (i-1) * 2, (i-1) * 2) = Eigen::MatrixXf::Identity((i-1) * 2, (i-1) * 2);
    a_nt.block(0, (i-1) * 2, (i-1) * 2, (i-1) * 2) = 2 * Eigen::MatrixXf::Identity((i-1) * 2, (i-1) * 2);

    time_res_operator.block(0, 0, (i-1) * 2, (i-1) * 2 * 2) = a_0;

    if(n%2 == 0){
        time_res_operator.block(time_res_operator.rows() - (i-1) * 2, time_res_operator.cols() - ((i-1) * 2 * 2), (i-1) * 2, (i-1) * 2 * 2) = a_nt;

        for(int j=1; j < n/2; ++j){

        time_res_operator.block(2*j*(i-1), 2*j*(i-1) + 6*(j-1), (i-1) * 2, (i-1) * 2 * 3) = a_m;

    }    
    }
    else{
        time_res_operator.block(time_res_operator.rows() - (i-1) * 2, time_res_operator.cols() - ((i-1) * 2), (i-1) * 2, (i-1) * 2) = Eigen::MatrixXf::Identity((i-1) * 2, (i-1) * 2);;
        
        for(int j=1; j < n/2 + 1; ++j){

        time_res_operator.block(2*j*(i-1), 2*j*(i-1) + 6*(j-1), (i-1) * 2, (i-1) * 2 * 3) = a_m;

    }
    }
    
    time_res_operator = (1.0/4) * time_res_operator;

    //std::cout << "\n time res operator: \n" << time_res_operator;

    d_time_restricted = time_res_operator * d;

    return d_time_restricted;
}
// Restriction operator ends here

// Restriction starts here
Eigen::VectorXf SpaceTimeRestriction(Eigen::VectorXf w_old, int n, int i){

    Eigen::VectorXf Time_w_res = Time_Restriction(w_old, n, i);
    //std::cout << "\n time_w_res size: \n" << Time_w_res.size();

    Eigen::VectorXf Space_w_res = Space_Restriction(Time_w_res, i);
    //std::cout << "\n w_res size: " << Space_w_res.size();

    return Space_w_res;
}

std::tuple<int, int, Eigen::VectorXf> Restriction(Eigen::VectorXf w_old, int n, int i){

    Eigen::VectorXf w_restricted = SpaceTimeRestriction(w_old, n, i);

    //std::cout << "w_restricted size: " << w_restricted.size() << std::endl;

    return std::make_tuple(n/2, i/2, w_restricted);
}
// Restriction ends here


// Prolongation starts here
// Space Prolongation starts here
Eigen::VectorXf Space_Prolongation(Eigen::VectorXf d, int i){

    Eigen::VectorXf d_space_prolongated;
    
    int n = d.size()/((i-1)*2);

    Eigen::MatrixXf space_prolongation = Eigen::MatrixXf::Zero( (2*(i-1) + 1) * 2 * n, d.size());

    Eigen::MatrixXf d_ = Eigen::MatrixXf::Zero(i-1, 2*(i-1) + 1);
    
    Eigen::MatrixXf d_m(1, 3);

    d_m << 1, 2, 1;
    
    for(int j=0; j < i-1; ++j){

        d_.block(j, 2*j, 1, 3) = d_m;
    }

    //std::cout << "\n d_: \n" << d_;

    Eigen::MatrixXf time_block = Eigen::MatrixXf::Zero((i-1) * 2, (2*(i-1) + 1)* 2);

    time_block.block(0, 0, i-1, 2*(i-1) + 1) = d_;
    time_block.block(i-1, 2*(i-1) + 1, i-1, 2*(i-1) + 1) = d_;

    Eigen::MatrixXf Prolongation_operator = time_block.transpose();

    //std::cout << "\n Prolongation_operator: \n" << Prolongation_operator;

    for(int j=0; j < n; ++j){

        space_prolongation.block((2*(i-1) + 1)* 2 * j, (i-1) * 2 * j, (2*(i-1) + 1)* 2, (i-1) * 2) = Prolongation_operator;
    }

    // std::cout << "\n space_prolongation: \n" << space_prolongation;

    d_space_prolongated = (1.0/4) * space_prolongation * d;

    return d_space_prolongated;
}
//Space Prolongation ends here

// Time prolongation starts here
Eigen::VectorXf Time_Prolongation(Eigen::VectorXf d, int i){

    Eigen::VectorXf d_time_prolongated;

    float n = d.size()/((i-1)*2);

    Eigen::MatrixXf time_prol_operator = Eigen::MatrixXf::Zero((n*2 - 1)*(i-1)*2, d.size());

    Eigen::MatrixXf time_block = Eigen::MatrixXf::Identity((i-1)*2, (i-1)*2);

    Eigen::MatrixXf time_block_even = Eigen::MatrixXf::Zero((i-1)*2, (i-1)*2*2);

    time_block_even.block(0, 0, (i-1)*2, (i-1)*2) = time_block;
    time_block_even.block(0, (i-1)*2, (i-1)*2, (i-1)*2) = time_block;

    for(int j=0; j < (n * 2 - 1); ++j){

        if(j%2 == 0){
            time_prol_operator.block((i-1)*2*j, (i-1)*j, (i-1)*2, (i-1)*2) = time_block;
        }
        
    }

    for(int j=0; j < (n * 2 - 1) - 1; ++j){

        if(j%2 == 0){
            time_prol_operator.block((i-1)*2*(j+1), (i-1)*j, (i-1)*2, (i-1)*2*2) = (1.0/2) * time_block_even;
        }
        
    }

    d_time_prolongated = (1.0/2) * time_prol_operator * d;

    return d_time_prolongated;
}
// Time prolongation ends here

// Combining time and space prolongation
Eigen::VectorXf SpaceTimeProlongation(Eigen::VectorXf w_old, int n, int i){

    Eigen::VectorXf Time_w_pro = Time_Prolongation(w_old, i);
    //std::cout << "\n time_w_prolongated size: \n" << Time_w_pro.size();

    Eigen::VectorXf Space_w_pro = Space_Prolongation(Time_w_pro, i);
    //std::cout << "\n w_prol size: " << Space_w_pro.size();

    return Space_w_pro;
}

std::tuple<int, int, Eigen::VectorXf> Prolongation(Eigen::VectorXf w_old, int n, int i){

    Eigen::VectorXf w_prol = SpaceTimeProlongation(w_old, n, i);

    //std::cout << "w_prol size: " << w_prol.size() << std::endl;

    return std::make_tuple(n*2, i*2, w_prol);
}
// Prolongation ends here

// Multilevel V cycle starts here
Eigen::VectorXf SpaceTimeMultigrid(Eigen::VectorXf w_old, Eigen::VectorXf f, int n, int i, float mu, int K, Eigen::VectorXf u){

    Eigen::VectorXf w_tilda, w_tilda_vector, w;
    Eigen::VectorXf error, d, residual;

    float delt = 1.0/n;
    float delx = 1.0/i;
    
    Eigen::MatrixXf A = BuildA(n, (i-1), delt, delx, mu, u);

    if (K==1){

        return SpaceTimeForwardBackwardBlockSORSmoother(A, w_old, f, n, (i-1), 20000, 1E-04);
        
    }
    
    w_tilda = SpaceTimeForwardBackwardBlockSORSmoother(A, w_old, f, n, (i-1), 4, 1E-08);    

    d = f - A * w_tilda;

    std::cout << "\n \n Restricting.... \n" << std::endl;

    auto d_tilda = Restriction(d, n, i);
    auto w_res = Restriction(w_tilda, n, i);

    float R_n = std::get<0>(d_tilda);
    float R_i = std::get<1>(d_tilda);
    delt = 1.0/R_n;
    delx = 1.0/R_i;

    Eigen::VectorXf u_res = Eigen::VectorXf::Zero((R_i-1)*(R_n+1));
    Eigen::VectorXf p_res = Eigen::VectorXf::Zero((R_i-1)*(R_n+1));
    Eigen::VectorXf y_res = Eigen::VectorXf::Zero((R_i-1)*(R_n+1));

    for(int j=0; j<R_n+1; j++){
        p_res.block(j*(R_i-1), 0, R_i-1, 1) = std::get<2>(w_res).block((2*j+1)*(R_i-1), 0, R_i-1, 1);
        y_res.block(j*(R_i-1), 0, R_i-1, 1) = std::get<2>(w_res).block((2*j*(R_i-1)), 0, R_i-1, 1);
    }

    for(int j=0; j<p_res.size(); j++){
        u_res(j) = (-1/mu) * p_res(j) * y_res(j);
    }

    Eigen::MatrixXf A_2h = BuildA(R_n, (R_i-1), delt, delx, mu, u_res);

    error = SpaceTimeMultigrid(Eigen::VectorXf::Zero(std::get<2>(d_tilda).size()), std::get<2>(d_tilda), R_n, R_i, mu, K-1, u_res);

    std::cout << "\n \n Prolongating.... \n" << std::endl;

    auto d_tilda_tilda_vector = Prolongation(error, std::get<0>(d_tilda), std::get<1>(d_tilda));

    w_tilda_vector =  w_tilda + std::get<2>(d_tilda_tilda_vector); 

    n = std::get<0>(d_tilda_tilda_vector);
    i = std::get<1>(d_tilda_tilda_vector);

    w = SpaceTimeForwardBackwardBlockSORSmoother(A, w_old, f, n, (i-1), 4, 1E-08); 

    return w;
}
// // Multilevel V cycle ends here

// Multiple Multilevel V cycle starts here
Eigen::VectorXf Multiple_SpaceTimeMultigrid(Eigen::VectorXf w_old, Eigen::VectorXf f, int n, int i, float mu, int K, Eigen::VectorXf u_old){

    std::ofstream out_file;

    remove("w_final.txt");

    Eigen::VectorXf w_tilda, w_tilda_vector, w;
    Eigen::VectorXf error, d, residual;

    float residual_norm = 1;
    float residual_norm_2 = 1;
    float tol = 0.01;

    float delt = 1.0/n;
    float delx = 1.0/i;
    int j = 0;

    Eigen::MatrixXf A = BuildA(n, (i-1), delt, delx, mu, u_old);

    if (K==1){
        
        return SpaceTimeForwardBackwardBlockSORSmoother(A, w_old, f, n, (i-1), 20000, 1E-04);
        
    }
    
    while((j < 4) && (residual_norm > tol)){
        
        std::cout << "\n Multigrid cycle: "  << j+1 << std::endl;
        std::cout << "\n Residual: "  << residual_norm << std::endl; 

        w_tilda = SpaceTimeForwardBackwardBlockSORSmoother(A, w_old, f, n, (i-1), 4, 1E-08);   

        d = f - A * w_tilda;

        std::cout << "\n \n Restricting.... \n" << std::endl;

        auto d_tilda = Restriction(d, n, i);
        auto w_res = Restriction(w_tilda, n, i);
        
        float R_n = std::get<0>(d_tilda);
        float R_i = std::get<1>(d_tilda);
        delt = 1.0/R_n;
        delx = 1.0/R_i;

        Eigen::VectorXf u_res = Eigen::VectorXf::Zero((R_i-1)*(R_n+1));
        Eigen::VectorXf p_res = Eigen::VectorXf::Zero((R_i-1)*(R_n+1));
        Eigen::VectorXf y_res = Eigen::VectorXf::Zero((R_i-1)*(R_n+1));

        for(int j=0; j<R_n+1; j++){
            p_res.block(j*(R_i-1), 0, R_i-1, 1) = std::get<2>(w_res).block((2*j+1)*(R_i-1), 0, R_i-1, 1);
            y_res.block(j*(R_i-1), 0, R_i-1, 1) = std::get<2>(w_res).block((2*j*(R_i-1)), 0, R_i-1, 1);
        }

        for(int j=0; j<p_res.size(); j++){
            u_res(j) = (-1/mu) * p_res(j) * y_res(j);
        }

        Eigen::MatrixXf A_2h = BuildA(R_n, (R_i-1), delt, delx, mu, u_res);

        error = SpaceTimeMultigrid(Eigen::VectorXf::Zero(std::get<2>(d_tilda).size()), std::get<2>(d_tilda), R_n, R_i, mu, K-1, u_res);

        std::cout << "\n \n Prolongating.... \n" << std::endl;

        auto d_tilda_tilda_vector = Prolongation(error, std::get<0>(d_tilda), std::get<1>(d_tilda));
        
        w_tilda_vector =  w_tilda + std::get<2>(d_tilda_tilda_vector);

        n = std::get<0>(d_tilda_tilda_vector);
        i = std::get<1>(d_tilda_tilda_vector);

        w = SpaceTimeForwardBackwardBlockSORSmoother(A, w_tilda_vector, f, n, (i-1), 4, 1E-08);

        residual_norm = (w - w_old).norm();

        w_old = w;
        
        j++;
    }

    std::cout << "\n Residual: "  << residual_norm << std::endl;   
    std::cout << "\n Multigrid V cycle Converged in: "  << j+1 << " Cycles " << std::endl; 
    
    out_file.open("multigrid_output.txt", std::ios_base::app); 
    out_file << "\n Multigrid V cycle Converged in: "  << j+1 << " Cycles "; 
    out_file.close();

    WriteFile(w, "w_final");

    return w;
}
// Multiple Multilevel V cycle ends here

// Collective smoothing multigrid
Eigen::VectorXf CSMG(Eigen::VectorXf w_old, Eigen::VectorXf f, int n, int i, float mu, int K, Eigen::VectorXf u_old){

    std::ofstream out_file;

    remove("w_final.txt");

    Eigen::VectorXf w_tilda, w_tilda_vector, w_vector, w;
    Eigen::VectorXf error, d, residual;

    float residual_norm = 1;
    float residual_norm_2 = 1;
    float tol = 0.1;

    float delt = 1.0/n;
    float delx = 1.0/i;
    int j = 0;

    Eigen::MatrixXf A = BuildA(n, (i-1), delt, delx, mu, u_old);

    if (K==1){
        
        return SpaceTimeForwardBackwardBlockSORSmoother(A, w_old, f, n, (i-1), 4, 1E-02);
        
    }
    
    while((residual_norm > tol)){
        
        std::cout << "\n Multigrid cycle: "  << j+1 << std::endl;
        std::cout << "\n Residual: "  << residual_norm << std::endl; 

        w_tilda = SpaceTimeForwardBackwardBlockSORSmoother(A, w_old, f, n, (i-1), 4, 1E-08);   

        // compute the residual
        d = f - A * w_tilda;

        std::cout << "\n \n Restricting.... \n" << std::endl;

        // restrict the residual
        auto d_tilda = Restriction(d, n, i);

        // restrict the solution
        auto w_res = Restriction(w_tilda, n, i);
        
        float R_n = std::get<0>(d_tilda);
        float R_i = std::get<1>(d_tilda);
        delt = 1.0/R_n;
        delx = 1.0/R_i;

        Eigen::VectorXf u_res = Eigen::VectorXf::Zero((R_i-1)*(R_n+1));
        Eigen::VectorXf p_res = Eigen::VectorXf::Zero((R_i-1)*(R_n+1));
        Eigen::VectorXf y_res = Eigen::VectorXf::Zero((R_i-1)*(R_n+1));

        for(int j=0; j<R_n+1; j++){
            p_res.block(j*(R_i-1), 0, R_i-1, 1) = std::get<2>(w_res).block((2*j+1)*(R_i-1), 0, R_i-1, 1);
            y_res.block(j*(R_i-1), 0, R_i-1, 1) = std::get<2>(w_res).block((2*j*(R_i-1)), 0, R_i-1, 1);
        }

        for(int j=0; j<p_res.size(); j++){
            u_res(j) = (-1/mu) * p_res(j) * y_res(j);
        }
 
        // compute the fine-to-coarse residual correction
        Eigen::MatrixXf A_2h = BuildA(R_n, (R_i-1), delt, delx, mu, u_res);
        auto f_res_old = Restriction(f, n, i);
        auto AW_res = Restriction(A * w_tilda, n, i);

        Eigen::VectorXf tau = (A_2h * std::get<2>(w_res)) - (std::get<2>(AW_res));
        Eigen::VectorXf f_res = (std::get<2>(f_res_old)) + tau;


        // CSMG to the coarse grid problem
        w_tilda_vector = CSMG(std::get<2>(w_res), f_res, R_n, R_i, mu, K-1, u_res);


        // error
        error = w_tilda_vector - std::get<2>(w_res);

        std::cout << "\n \n Prolongating.... \n" << std::endl;


        // prolongate the error
        auto error_pro = Prolongation(error, R_n, R_i);
        

        // coarse-grid corection
        w_vector =  w_tilda + std::get<2>(error_pro);


        // Post-smoothing
        w = SpaceTimeForwardBackwardBlockSORSmoother(A, w_vector, f, n, (i-1), 4, 1E-08);

        residual_norm = (w - w_old).norm();

        w_old = w;
        
        j++;
    }

    std::cout << "\n Residual: "  << residual_norm << std::endl;   
    std::cout << "\n Multigrid V cycle Converged in: "  << j << " Cycles " << std::endl; 
    
    out_file.open("multigrid_output.txt", std::ios_base::app); 
    out_file << "\n Multigrid V cycle Converged in: "  << j << " Cycles "; 
    out_file.close();

    WriteFile(w, "w_final");

    return w;
}
// Collective smoothing multigrid

// Single grid cycle
Eigen::VectorXf SpaceTimeSinglegrid(Eigen::VectorXf w_old, Eigen::VectorXf f, int n, int i, float mu, Eigen::VectorXf u_old){

    Eigen::VectorXf w, delta_w;
    
    float delt = 1.0/n;
    float delx = 1.0/i;

    float residual_norm = 1;
    float tol = 1e-03;

    int j=0;

    while((residual_norm >= tol)){
        
        std::cout << "Projection cycle: " << j+1 << std::endl;
        std::cout << "\n Projection cycle Residual: "  << residual_norm << std::endl;

        Eigen::MatrixXf A = BuildA(n, (i-1), delt, delx, mu, u_old);

        w = SpaceTimeForwardBackwardBlockSORSmoother(A, w_old, f, n, (i-1), 20000, 1E-03); 
        remove("w_internal.txt");
        WriteFile(w, "w_internal");

        Eigen::VectorXf u = Eigen::VectorXf::Zero((i-1)*(n+1));
        Eigen::VectorXf p = Eigen::VectorXf::Zero((i-1)*(n+1));
        Eigen::VectorXf y = Eigen::VectorXf::Zero((i-1)*(n+1));

        remove("p_internal.txt");
        remove("y_internal.txt");
        for(int j=0; j<n+1; j++){
            p.block(j*(i-1), 0, i-1, 1) = w.block((2*j+1)*(i-1), 0, i-1, 1);
            y.block(j*(i-1), 0, i-1, 1) = w.block((2*j*(i-1)), 0, i-1, 1);
        }
        WriteFile(p, "p_internal");
        WriteFile(y, "y_internal");

        remove("u.txt");
        for(int j=0; j<p.size() ; j++){
            u(j) = (-1/mu) * p(j) * y(j);
        }
        WriteFile(u, "u");

        residual_norm = (u - u_old).norm();

        w_old = w;
        u_old = u;
        j++;
}

std::cout << "\n Projection cycle Residual: "  << residual_norm << std::endl;   
std::cout << "\n Projection cycle Converged in: "  << j+1 << " Cycles" << std::endl; 
return w;

}
// Single grid cycle

// Multi grid cycle
Eigen::VectorXf NewtonProjectionMultigrid(Eigen::VectorXf w_old, Eigen::VectorXf f, int n, int i, float mu, Eigen::VectorXf u_old){

    Eigen::VectorXf w, delta_w;
    
    float delt = 1.0/n;
    float delx = 1.0/i;

    float residual_norm = 1;
    float tol = 0.1;

    int j=0;

    while((residual_norm >= tol)){
        
        std::cout << "Projection cycle: " << j+1 << std::endl;
        std::cout << "\n Projection cycle Residual: "  << residual_norm << std::endl;

        remove("w_internal.txt");
        //w = Multiple_SpaceTimeMultigrid(w_old, f, n, i, mu, 2, u_old);
        w = CSMG(w_old, f, n, i, mu, 2, u_old);
        WriteFile(w, "w_internal");

        Eigen::VectorXf u = Eigen::VectorXf::Zero((i-1)*(n+1));
        Eigen::VectorXf p = Eigen::VectorXf::Zero((i-1)*(n+1));
        Eigen::VectorXf y = Eigen::VectorXf::Zero((i-1)*(n+1));

        remove("p_internal.txt");
        remove("y_internal.txt");
        for(int j=0; j<n+1; j++){
            p.block(j*(i-1), 0, i-1, 1) = w.block((2*j+1)*(i-1), 0, i-1, 1);
            y.block(j*(i-1), 0, i-1, 1) = w.block((2*j*(i-1)), 0, i-1, 1);
        }
        WriteFile(p, "p_internal");
        WriteFile(y, "y_internal");

        remove("u.txt");
        for(int j=0; j<p.size() ; j++){
            u(j) = (-1/mu) * p(j) * y(j);
        }
        WriteFile(u, "u");

        residual_norm = (u - u_old).norm();

        w_old = w;
        u_old = u;
        j++;
}

std::cout << "\n Projection cycle Residual: "  << residual_norm << std::endl;   
std::cout << "\n Projection cycle Converged in: "  << j+1 << " Cycles" << std::endl; 
return w;

}
// Multi grid cycle

int main(){

    remove("w_final.txt");
    remove("u.txt");
    remove("Iterations.txt");
    remove("Residual.txt");

    int n = 128;
    int i = 64;
    
    float mu = 1;

    // Box bounded values
    float u_a = -0.1;
    float u_b = 0.1;

    float delt = 1.0/n;
    float delx = 1.0/i;

    std::cout << "\n delt: " << delt << std::endl;
    std::cout << "\n delx: " << delx << std::endl;

    // reading vectors Y, P
    Eigen::VectorXf Y = CSVtoVector("y.txt", i-1, n+1);
    Eigen::VectorXf P = CSVtoVector("p.txt", i-1, n+1);

    // reading vectors Z, F for unconstrained case
    Eigen::VectorXf z_unconst = CSVtoVector("z_unconst.txt", i-1, n+1);
    Eigen::VectorXf f_unconst = CSVtoVector("f_unconst.txt", i-1, n+1); 

    // Initializing w vector
    Eigen::VectorXf w_old = Eigen::VectorXf::Random(2*(i-1)*(n+1));

    // Initializing p and f vectors
    Eigen::VectorXf p = Eigen::VectorXf::Zero((i-1)*(n+1));
    Eigen::VectorXf y = Eigen::VectorXf::Zero((i-1)*(n+1));
    Eigen::VectorXf f = Eigen::VectorXf::Zero(2*(i-1)*(n+1)); 
        
    // applying boundary conditions and setting right hand side vector
    remove("f.txt");
    f = BoundaryConditions(f, Y, P, z_unconst, f_unconst, i-1, n);
    WriteFile(f, "f");

    remove("p_internal.txt");
    remove("y_internal.txt");
    Eigen::VectorXf u_old = Eigen::VectorXf::Zero((i-1)*(n+1));
    for(int j=0; j<n+1; j++){
        p.block(j*(i-1), 0, i-1, 1) = w_old.block((2*j+1)*(i-1), 0, i-1, 1);
        y.block(j*(i-1), 0, i-1, 1) = w_old.block((2*j*(i-1)), 0, i-1, 1);
    }
    WriteFile(p, "p_internal");
    WriteFile(y, "y_internal");

    remove("u.txt");
    for(int j=0; j<p.size(); j++){
        u_old(j) = (-1/mu) * p(j) * y(j);
    }
    WriteFile(u_old, "u");

    const clock_t begin_time = clock();
    std::ofstream output_file;
    remove("w.txt");
    Eigen::VectorXf w = NewtonProjectionMultigrid(w_old, f, n, i, mu, u_old);
    WriteFile(w, "w");
    output_file.open("Execution_time.txt", std::ios_base::app);
    output_file << "\n Execution time: "  << float( clock () - begin_time ) / CLOCKS_PER_SEC; 
    output_file.close();

    return 0;
}  
