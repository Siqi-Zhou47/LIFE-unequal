/*
    LIFE: Lattice boltzmann-Immersed boundary-Finite Element
    Copyright (C) 2019 Joseph O'Connor

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

// Includes
#include "../inc/FEMBody.h"
#include "../inc/Grid.h"
#include "../inc/Objects.h"
#include "../inc/Utils.h"
#include <sys/stat.h>
#include <sys/types.h>


// Dynamic FEM routine
void FEMBodyClass::dynamicFEM() {

	// Reset back to start of time step
	U = U_n;
	Udot = Udot_n;
	Udotdot = Udotdot_n;

	// Update FEM elements
	updateFEMValues();

	// Construct load vector as invariant during Newton-Raphson iterations
	constructRVector();

	// While loop parameters
	double TOL = 1e-10;
	double MAXIT = 20;

	// Set while counter to zero
	itNR = 0;

	// While loop for FEM solver
	do {

		// Solve and iterate over the system
		newtonRaphsonDynamic();

		// Check residual
		resNR = checkNRConvergence();

		// Increment counter
		itNR++;

	} while (resNR > TOL && itNR < MAXIT);

	// Compute new velocities and accelerations
	finishNewmark();

	// Update IBM nodes
	updateIBMValues();

	// Get subiteration residual for this body (will be summed later on)
	subResidual();
}

// Newton raphson iterator
void FEMBodyClass::newtonRaphsonDynamic() {

	// Build global matrices
	buildGlobalMatrices();

	// Apply Newmark scheme (using Newmark coefficients)
	setNewmark();

	// Solve linear system using LAPACK library
	delU = Utils::solveLAPACK(K, F, bcDOFs);

	// Add deltaU to U
	U = U + delU;

	// Update FEM positions
	updateFEMValues();
}

// Build global matrices
void FEMBodyClass::buildGlobalMatrices() {

	// Set matrices to zero
	fill(M.begin(), M.end(), 0.0);
	fill(K.begin(), K.end(), 0.0);
	fill(F.begin(), F.end(), 0.0);

	// Loop through and build global matrices
	for (size_t el = 0; el < element.size(); el++) {

		// Build force vector
		element[el].forceVector();

		// Build mass matrix
		element[el].massMatrix();

		// Build stiffness matrix
		element[el].stiffMatrix();
	}
}

// Set Newmark
void FEMBodyClass::setNewmark() {

	// Newmark-beta method for time integration
	double Dt = iPtr->oPtr->gPtr->Dt;
	double a0, a2, a3;
	a0 = 1.0 / (alpha * SQ(Dt));
	a2 = 1.0 / (alpha * Dt);
	a3 = 1.0 / (2.0 * alpha) - 1.0;

	// Calculate effective load vector
	F = R - F + Utils::MatMultiply(M, a0 * (U_n - U) + a2 * Udot + a3 * Udotdot);

	// Calculate effective stiffness matrix
	K = K + a0 * M;
}

// Finish Newmark
void FEMBodyClass::finishNewmark() {

	// Get timestep
	double Dt = iPtr->oPtr->gPtr->Dt;

	// Newmark coefficients
	double a6 = 1.0 / (alpha * SQ(Dt));
	double a7 = -1.0 / (alpha * Dt);
	double a8 = -(1.0 / (2.0 * alpha) - 1.0);
	double a9 = Dt * (1.0 - delta);
	double a10 = delta * Dt;

	// Update velocities and accelerations
	Udotdot = a6 * (U - U_n) + a7 * Udot_n + a8 * Udotdot_n;
	Udot = Udot_n + a9 * Udotdot_n + a10 * Udotdot;
}

// Update IBM markers
void FEMBodyClass::updateIBMValues() {

	// Declare local vectors
	array<double, elementDOFs> dashU;
	array<double, elementDOFs> dashUdot;
	array<double, dims> dashUShape;
	array<double, dims> dashUdotShape;
	array<array<double, dims>, dims> Tsub;

	// Loop through posMap points
	for (size_t i = 0; i < posMap.size(); i++) {

		// Get pointer to element and zeta value
		FEMElementClass *el = &(element[posMap[i].elID]);
		double zeta = posMap[i].zeta;

		// Disassemble positions, velocities and accelerations
		dashU = el->disassembleGlobalMat(U);
		dashUdot = el->disassembleGlobalMat(Udot);

		// Convert the displacement into local coordinates (shift angle as well)
		for (size_t n = 0; n < el->node.size(); n++) {
			dashU[n*nodeDOFs+eX] += el->node[n]->pos0[eX] - el->node[0]->pos[eX];
			dashU[n*nodeDOFs+eY] += el->node[n]->pos0[eY] - el->node[0]->pos[eY];
			dashU[n*nodeDOFs+(nodeDOFs-1)] += el->node[n]->angle0 - el->angle;
			dashU[n*nodeDOFs+(nodeDOFs-1)] = Utils::shiftAngle(dashU[n*nodeDOFs+(nodeDOFs-1)]);
		}

		// Get element values in local coordinates
		dashU = el->T * dashU;
		dashUdot = el->T * dashUdot;

		// Multiply by shape functions
		dashUShape = el->shapeFuns(dashU, zeta);
		dashUdotShape = el->shapeFuns(dashUdot, zeta);

		// Get subset of transformation matrix
		Tsub = {{{el->T[0][0], el->T[0][1]},
				 {el->T[1][0], el->T[1][1]}}};

		// Shift back to global coordinates
		dashUShape = Utils::Transpose(Tsub) * dashUShape;
		dashUdotShape = Utils::Transpose(Tsub) * dashUdotShape;

		// Set the IBM nodes
		iPtr->node[i]->pos = el->node[0]->pos + dashUShape;
		iPtr->node[i]->vel = dashUdotShape;
	}
}

// void FEMBodyClass::updateIBMValues() {
//     // If there is a posMap (mapping IBM markers to FEM elements + local zeta),
//     // use it. Otherwise fall back to simple along-beam linear interpolation.
//     size_t nIBM = iPtr->node.size();
//     if (nIBM == 0) return;

//     // Case A: posMap exists and matches IBM marker count
//     if (posMap.size() == nIBM) {
//         for (size_t i = 0; i < nIBM; ++i) {
//             const auto &pm = posMap[i];
//             // guard: ensure elID valid
//             if (pm.elID < 0 || pm.elID >= static_cast<int>(element.size())) {
//                 // fallback: put marker at base node
//                 iPtr->node[i]->pos = element.front().node[0]->pos;
//                 iPtr->node[i]->vel = {0.0, 0.0};
//                 continue;
//             }
//             FEMElementClass *el = &(element[pm.elID]);
//             double zeta = pm.zeta; // assumed in [-1,1]

//             // convert zeta (-1..1) to param s (0..1) for linear interpolation
//             double s = 0.5 * (zeta + 1.0);

//             // linear interpolation between el->node[0]->pos and el->node[1]->pos
//             array<double, dims> p0 = el->node[0]->pos;
//             array<double, dims> p1 = el->node[1]->pos;
//             array<double, dims> pos;
//             for (int d = 0; d < dims; ++d) pos[d] = (1.0 - s) * p0[d] + s * p1[d];

//             iPtr->node[i]->pos = pos;
//             // static case: zero velocity
//             iPtr->node[i]->vel = {0.0, 0.0};
//         }
//         return;
//     }

    // // Case B: no posMap or mismatch -> evenly distribute along whole filament (by arc-length between FEM nodes)
    // // We'll map each IBM marker index i to a normalized coordinate s in [0,1] along the filament,
    // // then find which FEM segment it belongs to and linearly interpolate there.
    // // Precompute cumulative element lengths
    // size_t nFEMNodes = node.size();
    // if (nFEMNodes < 2) {
    //     // degenerate: just set all markers to first node
    //     for (size_t i = 0; i < nIBM; ++i) {
    //         iPtr->node[i]->pos = node.front().pos;
    //         iPtr->node[i]->vel = {0.0, 0.0};
    //     }
    //     return;
    // }

//     // cumulative lengths along FEM nodes
//     vector<double> cumL(nFEMNodes, 0.0);
//     for (size_t k = 1; k < nFEMNodes; ++k) {
//         array<double, dims> v = { node[k].pos[0] - node[k-1].pos[0], node[k].pos[1] - node[k-1].pos[1] };
//         double Lk = sqrt(v[0]*v[0] + v[1]*v[1]);
//         cumL[k] = cumL[k-1] + Lk;
//     }
//     double totalL = cumL.back();
//     if (totalL <= 0.0) totalL = 1e-12; // avoid div0

//     for (size_t i = 0; i < nIBM; ++i) {
//         double s_global = (double)i / double(nIBM - 1); // 0..1
//         double targL = s_global * totalL;

//         // find element index idx such that cumL[idx] <= targL <= cumL[idx+1]
//         size_t idx = 0;
//         while (idx + 1 < nFEMNodes && cumL[idx+1] < targL) ++idx;
//         if (idx + 1 >= nFEMNodes) idx = nFEMNodes - 2;

//         double segL = cumL[idx+1] - cumL[idx];
//         double localS = (segL > 0.0) ? (targL - cumL[idx]) / segL : 0.0;
//         // linear interpolation between node[idx] and node[idx+1]
//         array<double, dims> p0 = node[idx].pos;
//         array<double, dims> p1 = node[idx+1].pos;
//         array<double, dims> pos;
//         for (int d = 0; d < dims; ++d) pos[d] = (1.0 - localS) * p0[d] + localS * p1[d];

//         iPtr->node[i]->pos = pos;
//         iPtr->node[i]->vel = {0.0, 0.0};
//     }
// }



// Update FEM values
void FEMBodyClass::updateFEMValues() {

	// Set the new positions
	for (size_t n = 0; n < node.size(); n++) {

		// Positions
		for (int d = 0; d < dims; d++)
			node[n].pos[d] = node[n].pos0[d] + U[node[n].DOFs[d]];

		// Set angle
		node[n].angle = node[n].angle0 + U[node[n].DOFs[dims]];
	}

	// Set the new angles and lengths of the elements
	array<double, dims> elVector;
	for (size_t el = 0; el < element.size(); el++) {

		// Set the new angles and lengths of the elements
		elVector = element[el].node[1]->pos - element[el].node[0]->pos;
		element[el].angle = atan2(elVector[eY], elVector[eX]);
		element[el].L = sqrt(elVector * elVector);

		// Set new transformation matrix
		element[el].setElementTransform();
	}
}

// Get load vector
void FEMBodyClass::constructRVector() {

	// Set R to zero
	fill(R.begin(), R.end(), 0.0);

	// Loop through all FEM elements
	for (size_t el = 0; el < element.size(); el++)
		element[el].loadVector();
}

// // Check convergence of Newton Raphson iterator
// inline double FEMBodyClass::checkNRConvergence () {

// 	// Get the norm of delU
// 	return sqrt(delU * delU) / (ref_L * sqrt(static_cast<double>(delU.size())));
// }

// void FEMBodyClass::constructRVector() {

// 	// Set R to zero
// 	fill(R.begin(), R.end(), 0.0);
// 	int dofPerNode = nodeDOFs;             // 一般为 3: (x, y, theta)
// 	 // Loop through all FEM elements (每个单元施加均布载荷)
//     if (qUniform != 0.0) {
//         double q = qUniform;                   // 均布载荷 [N/m]
//         double Le = L0 / (node.size() - 1);    // 单元长度
       

//         cout << "DEBUG: constructRVector, qUniform=" << qUniform
//              << " , totalForce=" << qUniform * L0 << " N" << endl;

//         for (size_t e = 0; e < node.size() - 1; ++e) {
//             // === 梁单元等效节点载荷 (Euler-Bernoulli) ===
//             double fe_yL     = q * Le / 2.0;
//             double fe_thetaL = q * Le * Le / 12.0;
//             double fe_yR     = q * Le / 2.0;
//             double fe_thetaR = -q * Le * Le / 12.0;

//             int leftBase  = e * dofPerNode;
//             int rightBase = (e + 1) * dofPerNode;

//             R[leftBase + 0]  += fe_yL;       // 左 y
//             R[leftBase + 2]  += fe_thetaL;   // 左 θ
//             R[rightBase + 0] += fe_yR;       // 右 y
//             R[rightBase + 2] += fe_thetaR;   // 右 θ
//         }
// 		// for (size_t n = 0; n < node.size(); ++n) {
//         //     int nodeBase = n * dofPerNode;
            
//         //     if (n == 0) {
//         //         // 底端节点：一半载荷
//         //         R[nodeBase + 0] += q * Le / 2.0;
//         //     }
//         //     else if (n == node.size() - 1) {
//         //         // 顶端节点：一半载荷
//         //         R[nodeBase + 0] += q * Le / 2.0;
//         //     }
//         //     else {
//         //         // 中间节点：全部载荷
//         //         R[nodeBase + 0] += q * Le;
//         //     }
//         // }
//     }
	
//     // 复制 K,R 到本地 (K 是 bodyDOFs x bodyDOFs 存扁平 vector, 行主序)
//     int n = bodyDOFs;
//     // 注意： 你可能需要先调用 buildGlobalMatrices() 来保证 K 已有值
//     vector<double> Kmat = K; // assume stored row-major
//     vector<double> Rvec = R;
//     vector<double> Utmp(n, 0.0);
//     // naive Gauss elimination (inplace on augmented matrix)
//     vector<double> A(n*(n+1),0.0);
//     for (int i=0;i<n;i++){
//       for (int j=0;j<n;j++) A[i*(n+1)+j] = Kmat[i*n + j];
//       A[i*(n+1)+n] = Rvec[i];
//     }
//     // forward elimination
//     for (int i=0;i<n;i++){
//       // find pivot
//       int piv = i;
//       double maxv = fabs(A[piv*(n+1)+i]);
//       for (int r=i+1;r<n;r++){
//         if (fabs(A[r*(n+1)+i]) > maxv){ maxv = fabs(A[r*(n+1)+i]); piv = r; }
//       }
//       if (maxv < 1e-16) continue; // singular
//       if (piv != i){
//         for (int c=i;c<=n;c++) swap(A[i*(n+1)+c], A[piv*(n+1)+c]);
//       }
//       double diag = A[i*(n+1)+i];
//       for (int r=i+1;r<n;r++){
//         double fac = A[r*(n+1)+i] / diag;
//         if (fac == 0.0) continue;
//         for (int c=i;c<=n;c++) A[r*(n+1)+c] -= fac * A[i*(n+1)+c];
//       }
//     }
//     // back substitution
//     for (int i=n-1;i>=0;i--){
//       double diag = A[i*(n+1)+i];
//       if (fabs(diag) < 1e-16){ Utmp[i] = 0.0; continue; } // singular -> skip
//       double s = A[i*(n+1)+n];
//       for (int c=i+1;c<n;c++) s -= A[i*(n+1)+c] * Utmp[c];
//       Utmp[i] = s / diag;
//     }

// 	// --- DIAGNOSTIC SNIPPET: paste after constructRVector() and before solve ---
// 	{
// 		std::ofstream diag("results/diag_constructRVector.txt", std::ios::app);
// 		int dofPerNode = nodeDOFs; // 3
// 		int nNodes = node.size();
// 		double L_sim = 0.1;         // total length you use in simulation
// 		double q_sim = 2.24e-3;   // uniform load you use (N/m)
// 		double E_sim = 4.0764;          // Young's modulus you use
// 		double D_sim = 0.01;

// 		// compute I
// 		double I = M_PI * pow(D_sim,4) / 64.0;

// 		// Sum of applied nodal loads (global)
// 		double sumRx=0.0, sumRy=0.0;
// 		for (int i=0;i<nNodes;i++){
// 			sumRx += R[i*dofPerNode + 0];
// 			sumRy += R[i*dofPerNode + 1];
// 		}

// 		// linear small-deflection estimate for tip x (if transverse load)
// 		// If your q is horizontal causing bending in global x, then treat q_t = q_sim
// 		double x_tip_lin = 0.0;
// 		// If bending is small and beam acts as cantilever under transverse load q (per unit length)
// 		// approximate tip displacement in transverse direction: delta = q * L^4 / (8 E I)
// 		// Use when problem corresponds to Euler-Bernoulli transverse load
// 		x_tip_lin = (q_sim * pow(L_sim,4)) / (8.0 * E_sim * I);

// 		// dimensionless stiffness K (paper): K = E I / (0.5 rho D L^3 U^2)
// 		double rho_f = 1000.0;
// 		double Uref = 0.02; double CD = 1.12; // only needed for Kpaper consistency
// 		double Kpaper = (E_sim * I) / (0.5 * rho_f * D_sim * pow(L_sim,3) * Uref * Uref);

// 		diag << "DIAG_CHECKS:\n";
// 		diag << "  L_sim="<<L_sim<<" D_sim="<<D_sim<<" E_sim="<<E_sim<<" I="<<I<<"\n";
// 		diag << "  q_sim="<<q_sim<<" sumRx="<<sumRx<<" sumRy="<<sumRy
// 			<<" expected q*L_x="<<q_sim*L_sim<<" expected q*L_y="<<0.0<<"\n";
// 		diag << "  linear small-deflection x_tip_lin="<<x_tip_lin
// 			<<" (non-dim x_tip/L ~ "<<(x_tip_lin/L_sim)<<")\n";
// 		diag << "  Kpaper="<<Kpaper<<" log10K="<<( (Kpaper>0)?log10(Kpaper):NAN )<<"\n";

// 		// print first few nodal rotations and positions
// 		diag << "  first nodes pos, theta(deg):\n";
// 		for (int i=0;i<std::min(8,(int)node.size()); ++i){
// 			double theta = U[i*dofPerNode + 2]; // if your U stores rotations
// 			diag << "    node "<<i<<" pos=("<<node[i].pos[0]<<","<<node[i].pos[1]<<")"
// 				<<" theta="<<theta*180.0/M_PI<<"\n";
// 		}
// 		diag << "---------------------------------------------\n";
//     	diag.close();
// 	}

// }


// void FEMBodyClass::constructRVector() {

//     fill(R.begin(), R.end(), 0.0);
//     int dofPerNode = nodeDOFs; // 一般 = 3
//     double Le = L0 / (node.size() - 1);
//     double q = qUniform; // 水平均布载荷[N/m]

//     cout << "DEBUG: constructRVector, qUniform=" << q
//          << ", totalForce=" << q * L0 << " N" << endl;

//     // 正确分布: 每个单元的等效节点载荷 (Euler-Bernoulli 梁)
//     for (size_t e = 0; e < node.size() - 1; ++e) {
//         int leftBase  = e * dofPerNode;
//         int rightBase = (e + 1) * dofPerNode;

//         // 力在 x 方向（水平），不产生纵向拉伸
//         double fe_xL = q * Le / 2.0;
//         double fe_xR = q * Le / 2.0;

//         // 注意：y方向不受力
//         R[leftBase + 0]  += fe_xL;
//         R[rightBase + 0] += fe_xR;
//     }

//     // 校验总力
//     double sumRx = 0.0;
//     for (size_t i = 0; i < R.size(); i += dofPerNode)
//         sumRx += R[i];
//     cout << "sumRx = " << sumRx << " expected q*L = " << q * L0 << endl;
// }


	// void FEMBodyClass::constructRVector() {
	// 	fill(R.begin(), R.end(), 0.0);
		
	// 	if (qUniform != 0.0) {
	// 		double q = qUniform;  // N/m，水平向右
			
	// 		for (size_t e = 0; e < element.size(); ++e) {
	// 			// 使用初始长度L0（不是当前长度L）
	// 			double Le = L0 / element.size();
				
	// 			// 全局X方向的均布载荷转换到当前单元局部坐标系
	// 			double theta = element[e].angle;
	// 			double q_local_x = q * cos(theta);  // 局部轴向分量
	// 			double q_local_y = - q * sin(theta);  // 局部横向分量
				
	// 			// 使用横向分量计算等效节点力（Euler-Bernoulli公式）
	// 			array<double, elementDOFs> R_local;
	// 			R_local[0] = q_local_x * Le / 2.0;        // 轴向力（会产生压缩/拉伸）
	// 			R_local[1] = q_local_y * Le / 2.0;        // 横向力
	// 			R_local[2] = q_local_y * Le * Le / 12.0;  // 弯矩
	// 			R_local[3] = q_local_x * Le / 2.0;
	// 			R_local[4] = q_local_y * Le / 2.0;
	// 			R_local[5] = -q_local_y * Le * Le / 12.0;
				
	// 			// 转换回全局坐标系
	// 			array<double, elementDOFs> R_global = Utils::Transpose(element[e].T) * R_local;
				
	// 			// 组装
	// 			element[e].assembleGlobalMat(R_global, R);
	// 		}
	// 	}
		
	// 	// 验证：总载荷应该等于 q * L0
	// 	double Rx_sum = 0, Ry_sum = 0, q = qUniform;
	// 	for (size_t i = 0; i < R.size(); i += nodeDOFs) {
	// 		Rx_sum += R[i];
	// 		Ry_sum += R[i+1];
	// 	}
	// 	if (fabs(Rx_sum - q * L0) > 1e-6 || fabs(Ry_sum) > 1e-6) {
	// 		cerr << "WARNING: Load sum incorrect! Rx=" << Rx_sum 
	// 			<< " (expect " << q*L0 << "), Ry=" << Ry_sum << endl;
	// 	}
	// }


// Replace existing constructRVector() with this function.
// Assumptions:
// - nodeDOFs = degrees per node (e.g. 3: x,y,theta)
// - element.size() == node.size()-1
// - element[e].L0 is initial element length (or global L0 / nElem)
// - element[e].angle0 is initial element angle (or use node pos0 to compute)
// - assembleGlobalMat(localVec, R) exists and adds to R
// void FEMBodyClass::constructRVector() {
//     fill(R.begin(), R.end(), 0.0);
//     if (qUniform == 0.0) return;

//     double q = qUniform; // global horizontal load [N/m], non-follower (fixed in space)
//     int dofPerNode = nodeDOFs;
//     size_t nElem = element.size();

//     // We'll compute per-element equivalent nodal forces using reference geometry (L0, angle0)
//     for (size_t e = 0; e < nElem; ++e) {
//         FEMElementClass &el = element[e];

//         // Use reference length (L0 per element). If element stores L0, use it; else compute L0 = L0_total / nElem
//         double Le = (el.L0 > 0.0) ? el.L0 : (L0 / double(nElem)); // safe fallback

//         // Compute element initial orientation using initial node positions if available
//         double x0_i = el.node[0]->pos0[eX];
//         double y0_i = el.node[0]->pos0[eY];
//         double x0_j = el.node[1]->pos0[eX];
//         double y0_j = el.node[1]->pos0[eY];
//         double theta0 = atan2(y0_j - y0_i, x0_j - x0_i); // initial tangent angle

//         // Global horizontal q vector is (q,0). We decompose into local components relative to initial element frame.
//         double c0 = cos(theta0), s0 = sin(theta0);
//         double q_local_axial = q * c0;     // axial component in initial local frame
//         double q_local_trans  = - q * s0;  // transverse component in initial local frame (sign chosen: +v is local up)

//         // Build local equivalent nodal vector (local DOFs: u1,v1,theta1, u2,v2,theta2)
//         std::array<double,6> Rloc = {0.0,0.0,0.0, 0.0,0.0,0.0};

//         // Axial contributions: equal half at each axial DOF
//         Rloc[0] += q_local_axial * Le / 2.0;
//         Rloc[3] += q_local_axial * Le / 2.0;

//         // Transverse contributions (standard Euler-Bernoulli consistent nodal loads for uniform transverse q_t)
//         // fy_left = q_t * L / 2; fy_right = q_t * L / 2
//         // M_left = q_t * L^2 / 12; M_right = - q_t * L^2 / 12
//         Rloc[1] += q_local_trans * Le / 2.0;     // v1
//         Rloc[4] += q_local_trans * Le / 2.0;     // v2
//         Rloc[2] +=  q_local_trans * Le * Le / 12.0;  // theta1
//         Rloc[5] += -q_local_trans * Le * Le / 12.0;  // theta2

//         // Transform local nodal loads (based on initial orientation) to global coordinates using el.T0 (initial transform)
//         // If you don't have T0 stored, form it here from theta0:
//         double T0[6][6] = {
//             { c0,  s0, 0,  0,  0, 0 },
//             {-s0,  c0, 0,  0,  0, 0 },
//             { 0,   0,  1,  0,  0, 0 },
//             { 0,   0,  0,  c0, s0, 0 },
//             { 0,   0,  0, -s0, c0, 0 },
//             { 0,   0,  0,  0,  0, 1 }
//         };

//         std::array<double,6> Rglob = {0,0,0,0,0,0};
//         for (int i=0;i<6;i++){
//             for (int j=0;j<6;j++){
//                 Rglob[i] += T0[j][i] * Rloc[j]; // Rglob = T0^T * Rloc
//             }
//         }

//         // Assemble into global R using element's assemble routine
//         el.assembleGlobalMat(Rglob, R);
//     }

//     // final sanity
//     double Rx_sum = 0.0, Ry_sum = 0.0;
//     for (size_t i = 0; i < R.size(); i += dofPerNode) {
//         Rx_sum += R[i];
//         Ry_sum += R[i+1];
//     }
//     cerr << "DEBUG: constructRVector (ref) sumRx = " << Rx_sum << " expected q*L0 = " << q*L0
//          << "  sumRy=" << Ry_sum << "\n";
// }



// Check convergence of Newton Raphson iterator
inline double FEMBodyClass::checkNRConvergence () {

	// Get the norm of delU
	return sqrt(delU * delU) / (ref_L * sqrt(static_cast<double>(delU.size())));
}

// Do a sum reduction to get the subiteration residual
void FEMBodyClass::subResidual() {

	// Get the residual from this time step and reassign old time step value
	R_km1 = R_k;
	R_k = U - U_km1;

	// Get residual for this body
	subRes = R_k * R_k;

	// Get numerator and denominator for calculating relaxation factor
	subNum = R_km1 * (R_k - R_km1);
	subDen = (R_k - R_km1) * (R_k - R_km1);
}

// predictor
void FEMBodyClass::predictor() {

	// Extrapolate the values
	if (iPtr->oPtr->gPtr->t > 2) {

		// Do 2rd order extrapolation
		U = 2.5 * U_n - 2.0 * U_nm1 + 0.5 * U_nm2;
	}
	else if (iPtr->oPtr->gPtr->t == 2) {

		// Do 1st order extrapolation
		U = 2.0 * U_n - U_nm1;
	}
	else if (iPtr->oPtr->gPtr->t == 1) {

		// Do zeroth order extrapolation
		U = U_n;
	}

	// Update FEM elements
	updateFEMValues();

	// Update the velocity
	finishNewmark();

	// Update IBM values
	updateIBMValues();

	// Set U_km1
	U_km1.swap(U);
}

// Get node mappings between FEM and IBM grids
void FEMBodyClass::computeNodeMapping(int nIBMNodes, int nFEMNodes) {

	// Size the map from FEM to IBM (1 for each IBM node)
	posMap.resize(nIBMNodes);

	// Set first node first
	posMap[0].elID = 0;
	posMap[0].zeta = -1.0;

	// Now loop through and set values
	for (int i = 1; i < nIBMNodes; i++) {

		// Check if remainder is zero
		if (fmod(i * ((nFEMNodes - 1.0) / (nIBMNodes - 1.0)), 1.0) == 0.0) {
			posMap[i].elID = static_cast<int>(std::floor(i * ((nFEMNodes - 1.0) / (nIBMNodes - 1.0)))) - 1;
			posMap[i].zeta = 1.0;
		}
		else {
			posMap[i].elID = static_cast<int>(std::floor(i * ((nFEMNodes - 1.0) / (nIBMNodes - 1.0))));
			posMap[i].zeta = fmod(i * ((nFEMNodes - 1.0) / (nIBMNodes - 1.0)), 1.0) * 2.0 - 1.0;
		}
	}

	// Loop through elements
	double node1, node2;
	for (size_t el = 0; el < element.size(); el++) {

		// Loop through all nodes and scale range to local coordinates for element
		for (int node = 0; node < nIBMNodes; node++) {
			node1 = -1 + 2.0 * ((static_cast<double>(node) - 0.5) * ((nFEMNodes - 1.0) / (nIBMNodes - 1.0)) - static_cast<double>(el)) / 1.0;
			node2 = -1 + 2.0 * ((static_cast<double>(node) + 0.5) * ((nFEMNodes - 1.0) / (nIBMNodes - 1.0)) - static_cast<double>(el)) / 1.0;

			// Check if any points lie within element coordinate range
			if ((node1 > -1.0 && node1 < 1.0) || (node2 > -1.0 && node2 < 1.0)) {

				// Sort out end nodes where one point lie outside element
				if (node1 < -1.0)
					node1 = -1.0;
				if (node2 > 1.0)
					node2 = 1.0;

				// Call constructor for chile IB point
				element[el].forceMap.emplace_back(node, node1, node2);
			}
		}
	}
}

// Reset the start of time step values
void FEMBodyClass::resetValues() {

	// Reset start of time step values
	U_nm2.swap(U_nm1);
	U_nm1.swap(U_n);
	U_n.swap(U);
	Udot_n.swap(Udot);
	Udotdot_n.swap(Udotdot);
}

// Set initial deflection using static FEM
void FEMBodyClass::setInitialDeflection() {

	// Write out
	cout << endl << "Starting static FEM for body " << iPtr->ID << "...";

	// Initial deflection
	double deflect = 0.0;

	// Set initial deflection
#ifdef INITIAL_DEFLECT
	deflect = -INITIAL_DEFLECT;
#endif

	// Get beam properties
	double E = element[0].E;
	double I = element[0].I;

	// Calculate linear force required to give specified deflection
	double linearForce = 3.0 * E * I * deflect / SQ(L0);

	// Set load steps and get max deltaForce
	int nSteps = 100;
	double deltaForceMax = linearForce / static_cast<double>(nSteps);

	// Get rotation matrix
	array<array<double, dims>, dims> T = Utils::getRotationMatrix(angle0);

	// Get initial force at tip
	array<double, dims> localForce = {0.0, deltaForceMax};
	array<double, dims> globalForce = T * localForce;

	// Set it to load vector
	R[bodyDOFs-nodeDOFs] = globalForce[0];
	R[bodyDOFs-nodeDOFs+1] = globalForce[1];

	// Declare feedback loop parameters
	array<double, dims> deflectVector;
	double gain = 1e-3;
	double TOL = 1e-10;
	int MAXIT = 10000;
	int it = 0;
	double error;

	// Start feedback loop
	do {

		// Static FEM calculation
		staticFEM();

		// Get deflection vector
		deflectVector =  {U[bodyDOFs-nodeDOFs], U[bodyDOFs-nodeDOFs+1]};
		deflectVector = Utils::Transpose(T) * deflectVector;

		// Get error
		error = (deflect - (deflectVector[1] / L0));

		// Force correction
		double forceCorrect = error * gain;

		// Check force correction is smaller than force interval
		if (fabs(forceCorrect) > fabs(deltaForceMax))
			forceCorrect = Utils::sgn(forceCorrect) * fabs(deltaForceMax);

		// Add to R vector
		localForce[1] += forceCorrect;
		globalForce = T * localForce;

		// Add to R vector
		R[bodyDOFs-nodeDOFs] = globalForce[0];
		R[bodyDOFs-nodeDOFs+1] = globalForce[1];

		// Increment counter
		it++;

		// Check if done max iterations
		if (it >= MAXIT)
			ERROR("max iterations hit! Try changing the gain...exiting");

	} while (fabs(error) > TOL);

	// Write out
	cout << "finished in " << it << " iterations";

	// Reset values to zero
	itNR = 0;
	resNR = 0.0;
	fill(R.begin(), R.end(), 0.0);
	fill(delU.begin(), delU.end(), 0.0);
}

// Newton raphson iterator
void FEMBodyClass::staticFEM() {

	// While loop parameters
	double TOL = 1e-10;
	double MAXIT = 20;

	// Set while counter to zero
	itNR = 0;

	// While loop for FEM solver
	do {

		// Solve and iterate over the system
		newtonRaphsonStatic();

		// Check residual
		resNR = checkNRConvergence();

		// Increment counter
		itNR++;

	} while (resNR > TOL && itNR < MAXIT);

	// Update IBM nodes
	updateIBMValues();
}


// void FEMBodyClass::staticFEM() {
//     cout << "\n=== [Static FEM Solver: Load-stepping Newton-Raphson] ===" << endl;

//     constructRVector();   // 构造全局载荷向量 R
//     const double TOL = 1e-8;       // 收敛容限
//     const int MAX_NR_IT = 200;     // 每个载荷步的最大NR迭代次数
//     const int N_LOAD_STEPS = 100;  // 总载荷分步数（越多越稳）
//     const double loadScale = 1.0 / N_LOAD_STEPS;

//     fill(U.begin(), U.end(), 0.0); // 初始位移全零
//     std::vector<double> Rstep(R.size(), 0.0); // 当前步的目标载荷
//     double totalResidual = 0.0;

//     namespace fs = boost::filesystem;
//     fs::create_directories("Staticresults");

// 	int base = 0 * nodeDOFs;  // 根节点的自由度索引起点

// 	// === 强制根节点固定 ===
// 	for (int d = 0; d < nodeDOFs; ++d) {
// 		U[base + d] = 0.0;
// 		F[base + d] = 0.0;
// 		R[base + d] = 0.0;
// 	}

//     // === Load stepping ===
//     for (int loadStep = 1; loadStep <= N_LOAD_STEPS; ++loadStep) {
//         cout << "\n[Load Step " << loadStep << "/" << N_LOAD_STEPS << "]" << endl;
//         double factor = loadStep * loadScale;

//         // 当前载荷比例
//         for (size_t i = 0; i < R.size(); ++i)
//             Rstep[i] = factor * R[i];

//         // Newton-Raphson 初始化
//         int itNR = 0;
//         double resNR = 1e10;
// 		for (int d = 0; d < nodeDOFs; ++d) {
// 			U[base + d] = 0.0;
// 			F[base + d] = 0.0;
// 			R[base + d] = 0.0;
// 		}

//         // === Newton-Raphson 循环 ===
//         while (resNR > TOL && itNR < MAX_NR_IT) {
//             itNR++;
			
// 			//   // 更新位移场
//             // U = U + delU;
//             // updateFEMValues();
// 			// //constructRVector(); 
//             // updateIBMValues();

// 			buildGlobalMatrices();
// 			// --- DIAGNOSTIC: 在 solve 之前检查 K, R 是否包含 NaN/Inf 和奇异值 ---
// 			bool bad = false;
// 			for (int i=0;i<bodyDOFs;i++){
// 				double rv = R[i];
// 				if (!isfinite(rv)) {
// 					cerr << "[ERROR] R["<<i<<"] is not finite: " << rv << "\n";
// 					bad = true;
// 					break;
// 				}
// 			}
// 			for (int i=0;i<bodyDOFs && !bad;i++){
// 				for (int j=0;j<bodyDOFs;j++){
// 					double kv = K[i*bodyDOFs + j];
// 					if (!isfinite(kv)) {
// 						cerr << "[ERROR] K["<<i<<","<<j<<"] is not finite: " << kv << "\n";
// 						bad = true;
// 						break;
// 					}
// 				}
// 			}
// 			if (bad) {
// 				// Dump some diagnostics and abort safely
// 				std::ofstream of("Staticresults/debug_nan_dump.txt", std::ios::app);
// 				of << "STEP " << itNR << " -> NaN/Inf detected\n";
// 				of << "First few R: ";
// 				for (int i=0;i<std::min(bodyDOFs,10);++i) of << R[i] << " ";
// 				of << "\nK diag (first 10): ";
// 				for (int i=0;i<std::min(bodyDOFs,10);++i) of << K[i*bodyDOFs + i] << " ";
// 				of << "\n";
// 				of.close();
// 				ERROR("NaN/Inf in K or R detected - aborting static solve (see debug_nan_dump.txt)");
// 			}

//             // 解增量方程： K * delU = Rstep - F
//             delU = Utils::solveLAPACK(K, Rstep - F, bcDOFs);
// 			// check delU for NaN/Inf
// 			for (size_t i=0;i<delU.size();++i){
// 				if (!isfinite(delU[i])) {
// 					cerr << "[ERROR] delU["<<i<<"] is not finite: "<< delU[i] << "\n";
// 					// dump K diag and small portion of K/R for debugging
// 					std::ofstream of("Staticresults/solver_bad_dump.txt", std::ios::app);
// 					of << "Iteration " << itNR << " delU contains non-finite\n";
// 					of << "First 10 R: ";
// 					for (int k=0;k<std::min((int)R.size(),10);++k) of << R[k] << " ";
// 					of << "\nK diag: ";
// 					for (int k=0;k<std::min(bodyDOFs,10);++k) of << K[k*bodyDOFs + k] << " ";
// 					of << "\n";
// 					of.close();
// 					ERROR("Solver returned non-finite delU - aborting");
// 				}
// 			}

//             // 更新位移场
//             U = U + delU;
//             updateFEMValues();
// 			//constructRVector(); 
//             updateIBMValues();

//             // 检查收敛
//             resNR = checkNRConvergence();
//             cout << "  Iteration " << itNR << ", residual = " << resNR << endl;

//             if (itNR % 50 == 0) {
//                 cout << "  [debug] mid-step residual check at iteration " << itNR << endl;
//             }

//             if (std::isnan(resNR) || resNR > 1e8) {
//                 cerr << "  WARNING: Divergence detected! Breaking..." << endl;
//                 break;
//             }
//         }

//         cout << "  --> Converged in " << itNR << " iterations (res=" << resNR << ")\n";

//         // === 写出VTK文件 ===
//         if (loadStep % 10 == 0 || loadStep == N_LOAD_STEPS) {
//             std::string vtkName = "Staticresults/static_loadstep_" + std::to_string(loadStep) + ".vtk";
//             std::ofstream vtk(vtkName);
//             vtk << "# vtk DataFile Version 3.0\n";
//             vtk << "Static load step " << loadStep << "\n";
//             vtk << "ASCII\n";
//             vtk << "DATASET POLYDATA\n";
//             vtk << "POINTS " << iPtr->node.size() << " float\n";
//             for (auto n : iPtr->node)
//                 vtk << n->pos[0] << " " << n->pos[1] << " 0.0\n";
//             vtk << "LINES " << (iPtr->node.size()-1) << " " << 3*(iPtr->node.size()-1) << "\n";
//             for (size_t i=0;i<iPtr->node.size()-1;i++)
//                 vtk << "2 " << i << " " << i+1 << "\n";
//             vtk.close();
//             cout << "  [VTK] Wrote " << vtkName << endl;
//         }

//         totalResidual += resNR;
//         if (resNR > 1.0) {
//             cout << "  [warn] residual still large at step " << loadStep << endl;
//         }
//     }

//     // Put this right after constructRVector();
// 	{
// 		std::cerr << "=== QUICK DIAG START ===\n";
// 		// sums
// 		double sumRx=0, sumRy=0;
// 		for (int i=0;i<(int)R.size(); i+=nodeDOFs) { sumRx += R[i+0]; sumRy += R[i+1]; }
// 		std::cerr << "sumRx = " << sumRx << "  sumRy = " << sumRy << "  expected q*L = " << qUniform*L0 << "\n";

// 		// print first few R
// 		std::cerr << "R first entries: ";
// 		for (int i=0;i<std::min((int)R.size(), 12); ++i) std::cerr << R[i] << " ";
// 		std::cerr << "\n";

// 		// K diag
// 		double avgKdiag=0;
// 		int cnt = std::min(10, bodyDOFs);
// 		for (int i=0;i<cnt;i++) avgKdiag += fabs(K[i*bodyDOFs + i]);
// 		avgKdiag /= std::max(1,cnt);
// 		std::cerr << "avgKdiag = " << avgKdiag << " ; first diag: ";
// 		for (int i=0;i<cnt;i++) std::cerr << K[i*bodyDOFs + i] << " ";
// 		std::cerr << "\n";

// 		std::cerr << "DIAG_CHECKS: bcDOFs (scalar) = " << bcDOFs << "\n";

// 		// material/geometric analytic check
// 		double Dsim = 0.01;
// 		double Lsim = L0; // ensure correct
// 		double Esim = 407.64;
// 		double Isim = M_PI * pow(Dsim,4)/64.0;
// 		double EI = Esim * Isim;
// 		double qsim = qUniform;
// 		double xtip_lin = qsim * pow(Lsim,4) / (8.0 * EI);
// 		double tipSlope_deg = (qsim * pow(Lsim,3) / (6.0 * EI)) * 180.0 / M_PI;
// 		std::cerr << "ANALYTIC: D="<<Dsim<<" L="<<Lsim<<" E="<<Esim<<" I="<<Isim<<" EI="<<EI
// 				<< " q="<<qsim<<" xtip_lin="<<xtip_lin<<" tipSlope(deg)="<<tipSlope_deg<<"\n";
// 		std::cerr << "=== QUICK DIAG END ===\n";
// 	}

// 				// === 完整的能量检查 ===
// 			double strain_energy = 0.0;
// 			double work_external = 0.0;
			
// 			// 计算应变能
// 			for (size_t e = 0; e < element.size(); e++) {
// 				double u = (element[e].L - element[e].L0);
// 				double theta1 = Utils::shiftAngle(element[e].node[0]->angle - element[e].angle);
// 				double theta2 = Utils::shiftAngle(element[e].node[1]->angle - element[e].angle);
				
// 				double E = element[e].E;
// 				double A = element[e].A;
// 				double I = element[e].I;
// 				double L0 = element[e].L0;
				
// 				// 轴向应变能
// 				double U_axial = 0.5 * (E * A / L0) * u * u;
				
// 				// 弯曲应变能
// 				double U_bending = (E * I / L0) * (2.0*theta1*theta1 + 2.0*theta2*theta2 + 2.0*theta1*theta2);
				
// 				strain_energy += U_axial + U_bending;
// 			}
			
// 			// 计算外力功
// 			for (size_t i = 0; i < U.size(); i++) {
// 				work_external += R[i] * U[i];
// 			}
			
// 			cout << "\n=== ENERGY CHECK ===" << endl;
// 			cout << "Strain energy: " << strain_energy << " J" << endl;
// 			cout << "External work: " << work_external << " J" << endl;
// 			cout << "Ratio (should be ~0.5 for linear elastic): " << strain_energy / work_external << endl;
			
// 			// 如果比值不是 0.5 左右，说明有问题
// 			if (fabs(strain_energy / work_external - 0.5) > 0.1) {
// 				cerr << "[WARNING] Energy imbalance detected!" << endl;
// 				cerr << "This suggests inconsistency between stiffness and force calculation!" << endl;
// 			}
			
// 			// === 手动计算理论值（精确积分）===
// 			double q = qUniform;
// 			double E_mat = element[0].E;
// 			double I_mat = element[0].I;
// 			double L_tot = L0;
			
// 			// 使用精确的大变形理论（椭圆积分方法）
// 			// 或者使用 Bisshopp-Drucker 大挠度悬臂梁表
			
// 			cout << "\n=== LARGE DEFLECTION THEORY ===" << endl;
// 			cout << "For comparison, need to use large deflection beam theory" << endl;
// 			cout << "(Bisshopp-Drucker tables or elliptic integrals)" << endl;
			
// 			// 无量纲参数
// 			double Lambda = q * pow(L_tot, 3) / (E_mat * I_mat);
// 			cout << "Load parameter Lambda = qL³/(EI) = " << Lambda << endl;
			
// 			// 从 Bisshopp-Drucker 表（如果 Lambda 很大，需要用大变形理论）
// 			// Lambda > 1 时，小变形理论误差 > 10%
			
// 			if (Lambda > 1.0) {
// 				cout << "[INFO] Lambda > 1, this is a LARGE deflection problem!" << endl;
// 				cout << "Small deflection theory error > 10%" << endl;
// 				cout << "Need to compare with Bisshopp-Drucker tables or exact solution" << endl;
// 			}
// 		cout << "\n=== Static analysis finished ===" << endl;
// 		cout << "Total accumulated residual = " << totalResidual << endl;
// 		auto &tipNode = node.back();
// 		double tipDx = tipNode.pos[0] - node.front().pos[0];
// 		double tipDy = tipNode.pos[1] - node.front().pos[1];
// 		cout << ">>> Tip displacement: Dx = " << tipDx << ", Dy = " << tipDy << endl;
// 		double angle_rad = atan2(tipDy, tipDx);
// 		double angle_deg = angle_rad * 180.0 / M_PI;
// 		cout << ">>> Tip angle = " << angle_deg << " deg" << endl;

// }

// void FEMBodyClass::staticFEM() {

//     cerr << "DEBUG IBMBody: running staticFEM() for standalone static deformation...\n";

//     const int N_LOAD_STEPS = 500;
//     const int MAX_ITER = 100;
//     const double TOL = 1e-8;

//     double totalAccumResidual = 0.0;

//     fill(U.begin(), U.end(), 0.0);

//     for (int step = 1; step <= N_LOAD_STEPS; ++step) {

//         double loadFactor = double(step) / N_LOAD_STEPS;
//         constructRVector_scaled(loadFactor);

//         double residual = 1.0;
//         int iter = 0;

//         do {
//             buildGlobalMatrices();
//             delU = Utils::solveLAPACK(K, R - F, bcDOFs);
//             U = U + delU;
//             updateFEMValues();
// 			updateIBMValues();
//             residual = checkNRConvergence();
//             iter++;
//         } while (residual > TOL && iter < MAX_ITER);

//         totalAccumResidual += residual;

//         if (step % 10 == 0) {
//             cerr << "  [Load Step " << step << "/" << N_LOAD_STEPS << "] "
//                  << "Converged in " << iter << " iterations (res=" << residual << ")\n";

//             namespace fs = boost::filesystem;
//             fs::create_directories("Staticresults");
//             std::ofstream vtk("Staticresults/static_loadstep_" + std::to_string(step) + ".vtk");
//             vtk << "# vtk DataFile Version 3.0\n";
//             vtk << "Static load step " << step << "\n";
//             vtk << "ASCII\n";
//             vtk << "DATASET POLYDATA\n";
//             vtk << "POINTS " << iPtr->node.size() << " float\n";
//             for (auto n : iPtr->node)
//                 vtk << n->pos[0] << " " << n->pos[1] << " 0.0\n";
//             vtk << "LINES " << (iPtr->node.size()-1) << " " << 3*(iPtr->node.size()-1) << "\n";
//             for (size_t i=0;i<iPtr->node.size()-1;i++)
//                 vtk << "2 " << i << " " << i+1 << "\n";
//             vtk.close();
//         }

//         if (residual > 1.0) {
//             cerr << "  WARNING: Divergence detected at step " << step << "!\n";
//             break;
//         }
//     }

//     cerr << "\n=== Static analysis finished ===\n";
//     cerr << "Total accumulated residual = " << totalAccumResidual << "\n";
// 		auto &tipNode = node.back();
// 	double tipDx = tipNode.pos[0] - node.front().pos[0];
// 	double tipDy = tipNode.pos[1] - node.front().pos[1];
// 	cout << ">>> Tip displacement: Dx = " << tipDx << ", Dy = " << tipDy << endl;
// 	double angle_rad = atan2(tipDy, tipDx);
// 	double angle_deg = angle_rad * 180.0 / M_PI;
// 	cout << ">>> Tip angle = " << angle_deg << " deg" << endl;
// }


// Newton raphson iterator
void FEMBodyClass::newtonRaphsonStatic() {

	// Build global matrices
	buildGlobalMatrices();

	// Solve linear system using LAPACK library
	delU = Utils::solveLAPACK(K, R - F, bcDOFs);

	// Add deltaU to U
	U = U + delU;

	// Update FEM positions
	updateFEMValues();
}

// Custom constructor for building corotational FEM body
FEMBodyClass::FEMBodyClass(IBMBodyClass *iBodyPtr, const array<double, dims> &pos, const array<double, dims> &geom, double angle, string nElementsStr, string BC, double rho, double E, double qUniform_param) {

	// Set pointer
	iPtr = iBodyPtr;
	
	// Set to initial value
	itNR = 0;
	resNR = 0.0;

	// Set sub residaul, numerator and denominator to initial value
	subRes = 0.0;
	subNum = 0.0;
	subDen = 0.0;
	this->qUniform = qUniform_param;

	// Get number of DOFs required for BC
	if (BC == "CLAMPED")
		bcDOFs = 3;
	else if (BC == "SUPPORTED")
		bcDOFs = 2;

	// Unpack geometry
	angle0 = angle;
	L0 = geom[0];

	// Get number of elements
	int numElements;
	if (nElementsStr == "CONFORMING")
		numElements = static_cast<int>(floor(L0 / iPtr->oPtr->gPtr->Dx));
	else
		numElements = static_cast<int>(stod(nElementsStr));

	// Get number of nodes
	int numNodes = numElements + 1;

	// Get number of DOFs in body
	bodyDOFs = numNodes * nodeDOFs;

	// Get critical time step for this body
	tCrit = (L0 / numElements) / sqrt(E / rho);

	// Get rotation matrix
	array<array<double, dims>, dims> T = Utils::getRotationMatrix(angle);

	// Position vector for marker
	array<double, dims> position;

	// Loop through and build nodes
	for (int i = 0; i < numNodes; i++) {

		// Get position of this marker
		position[eX] = i * L0 / numElements;
		position[eY] = 0.0;

		// Rotate
		position = T * position;

		// Add the start point
		position = pos + position;

		// Call node constructor
		node.emplace_back(i, position, angle);
	}
	
	// Loop through and build elements
	for (int i = 0; i < numElements; i++){
		element.emplace_back(this, i, geom, angle, L0 / numElements, rho, E);
		
	}
		
	// Get number of IBM and FEM nodes
	int nIBMNodes = static_cast<int>(iPtr->node.size());
	int nFEMNodes = static_cast<int>(node.size());

	// Compute IBM-FEM conforming parameters
	computeNodeMapping(nIBMNodes, nFEMNodes);

	// Size the matrices
	M.resize(bodyDOFs * bodyDOFs, 0.0);
	K.resize(bodyDOFs * bodyDOFs, 0.0);
	R.resize(bodyDOFs, 0.0);
	F.resize(bodyDOFs, 0.0);
	U.resize(bodyDOFs, 0.0);
	delU.resize(bodyDOFs, 0.0);
	Udot.resize(bodyDOFs, 0.0);
	Udotdot.resize(bodyDOFs, 0.0);
	U_n.resize(bodyDOFs, 0.0);
	U_nm1.resize(bodyDOFs, 0.0);
	U_nm2.resize(bodyDOFs, 0.0);
	Udot_n.resize(bodyDOFs, 0.0);
	Udotdot_n.resize(bodyDOFs, 0.0);
	U_km1.resize(bodyDOFs, 0.0);
	R_k.resize(bodyDOFs, 0.0);
	R_km1.resize(bodyDOFs, 0.0);
	computeNodeMapping(nIBMNodes, nFEMNodes);
}
