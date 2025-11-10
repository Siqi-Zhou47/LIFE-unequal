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
#include "../inc/FEMElement.h"
#include "../inc/FEMBody.h"
#include "../inc/Grid.h"
#include "../inc/Objects.h"
#include "../inc/Utils.h"


// Construct load vector
void FEMElementClass::loadVector() {

	// Declare vectors
	array<array<double, dims>, dims> Tsub;
	array<double, dims> F;
	array<double, elementDOFs> RGlobal;

	// Get weight vector
	array<double, dims> weight = {rho * A * gravityX, rho * A * gravityY};

	// Get force scaling parameter
	double forceScale = fPtr->iPtr->oPtr->gPtr->Dm / SQ(fPtr->iPtr->oPtr->gPtr->Dt);

	// Loop through all force mapping nodes
	for (size_t n = 0; n < forceMap.size(); n++) {

		// Get node and integration ranges
		IBMNodeClass *node = fPtr->iPtr->node[forceMap[n].nodeID];
		double a = forceMap[n].zeta1;
		double b = forceMap[n].zeta2;

		// Get subset of transpose matrix
		Tsub = {{{T[0][0], T[0][1]},
				 {T[1][0], T[1][1]}}};

		// Convert force to local coordinates
		F = Tsub * (((-node->epsilon * 1.0 * forceScale) * node->force) + weight);

		// Get the nodal values by integrating over range of IB point
		R[0] = F[0] * 0.5 * L * (0.5 * b - 0.5 * a + 0.25 * SQ(a) - 0.25 * SQ(b));
		R[1] = F[1] * 0.5 * L * (0.5 * b - 0.5 * a - SQ(a) * SQ(a) / 16.0 + SQ(b) * SQ(b) / 16.0 + 3.0 * SQ(a) / 8.0 - 3.0 * SQ(b) / 8.0);
		R[2] = F[1] * 0.5 * L * (L * (-SQ(a) * SQ(a) + SQ(b) * SQ(b)) / 32.0 - L * (-TH(a) + TH(b)) / 24.0 - L * (-SQ(a) + SQ(b)) / 16.0 + L * (b - a) / 8.0);
		R[3] = F[0] * 0.5 * L * (-0.25 * SQ(a) + 0.25 * SQ(b) + 0.5 * b - 0.5 * a);
		R[4] = F[1] * 0.5 * L * (0.5 * b - 0.5 * a + SQ(a) * SQ(a) / 16.0 - SQ(b) * SQ(b) / 16.0 - 3.0 * SQ(a) / 8.0 + 3.0 * SQ(b) / 8.0);
		R[5] = F[1] * 0.5 * L * (L * (-SQ(a) * SQ(a) + SQ(b) * SQ(b)) / 32.0 + L * (-TH(a) + TH(b)) / 24.0 - L * (-SQ(a) + SQ(b)) / 16.0 - L * (b - a) / 8.0);

		// Get element load vector
		RGlobal = Utils::Transpose(T) * R;

		// Assemble into global vector
		assembleGlobalMat(RGlobal, fPtr->R);
	}
}

// Construct internal force vector
void FEMElementClass::forceVector() {

	// Calculate the local displacements
	double u = (SQ(L) - SQ(L0)) / (L + L0);
	
	double theta1 = Utils::shiftAngle(node[0]->angle - angle);
	double theta2 = Utils::shiftAngle(node[1]->angle - angle);

	// Calculate internal forces for beam
	double F0 = (E * A / L0) * u;
	double M1 = (2 * E * I / L0) * (2.0 * theta1 + theta2);
	double M2 = (2 * E * I / L0) * (theta1 + 2.0 * theta2);

	// Set the internal local nodal forces
	F[0] = -F0;
	F[1] = (1.0 / L0) * (M1 + M2);
	F[2] = M1;
	F[3] = F0;
	F[4] = -(1.0 / L0) * (M1 + M2);
	F[5] = M2;

	// Get element internal forces
	array<double, elementDOFs> FGlobal = Utils::Transpose(T) * F;
	
	// Assemble into global vector
	assembleGlobalMat(FGlobal, fPtr->F);
}

	// void FEMElementClass::forceVector() {

	//     // --- 1) 更新当前几何（端点位置） ---
	//     double xi = node[0]->pos[0];
	//     double yi = node[0]->pos[1];
	//     double xj = node[1]->pos[0];
	//     double yj = node[1]->pos[1];

	//     L = sqrt((xj - xi)*(xj - xi) + (yj - yi)*(yj - yi));
	//     if (!isfinite(L) || L < 1e-12) {
	// 		cerr << "[WARN] element " << this->elemID << " tiny or NaN length L="<<L<<", using L0 fallback\n";
	// 		L = L0;
	// 	}
	//     c = (xj - xi) / L;
	//     s = (yj - yi) / L;
	// 	if (!isfinite(c) || !isfinite(s)) {
	// 		cerr << "[ERROR] non-finite direction cosines c="<<c<<" s="<<s<<" at elem="<<this->elemID<<"\n";
	// 		c = 1.0; s = 0.0; // fallback (aligned with x) to avoid NaN propagation
	// 	}
		
	// 	// ensure Flocal initialized
	// 	for (int i=0;i<6;i++) Flocal[i] = 0.0;
	//     // --- 2) 计算局部化位移/角度（基于你原实现的思想，但注意用当前几何） ---
	//     // u: axial relative extension (you had (L^2-L0^2)/(L+L0) earlier)
	//     // double u = (L*L - L0*L0) / (L + L0);
	// 	double log_strain = log(L / L0);
	// 	double F0 = E * A * log_strain;
	//     // theta offset: 使用节点相对单元 reference angle (原来 node[i]->angle - angle)
	//     double theta1 = Utils::shiftAngle(node[0]->angle - angle);
	//     double theta2 = Utils::shiftAngle(node[1]->angle - angle);

	//     // --- 3) 局部内力（轴力和弯矩）按原公式（保持一致） ---
	
	//     double M1 = (2.0 * E * I / L) * (2.0 * theta1 + theta2);
	//     double M2 = (2.0 * E * I / L) * (theta1 + 2.0 * theta2);

	//     // Flocal: 局部节点内力向量 (uL, vL, mL, uR, vR, mR)
	//     Flocal[0] = -F0;
	//     Flocal[1] = (1.0 / L) * (M1 + M2);
	//     Flocal[2] = M1;
	//     Flocal[3] = F0;
	//     Flocal[4] = -(1.0 / L) * (M1 + M2);
	//     Flocal[5] = M2;

	//     // --- 4) 将局部力旋转到全局并装配 ---
	//     // 转换矩阵 T (6x6)
	//     double Tmat[6][6] = {
	//         { c,  s, 0,  0,  0, 0 },
	//         { -s, c, 0,  0,  0, 0 },
	//         { 0,  0, 1,  0,  0, 0 },
	//         { 0,  0, 0,  c,  s, 0 },
	//         { 0,  0, 0, -s,  c, 0 },
	//         { 0,  0, 0,  0,  0, 1 }
	//     };

	//     // Flocal -> Fglobal: Fg = T^T * Flocal
	//     std::array<double,6> Fglobal;
	//     for (int i=0;i<6;i++){
	//         Fglobal[i] = 0.0;
	//         for (int j=0;j<6;j++){
	//             Fglobal[i] += Tmat[j][i] * Flocal[j];
	//         }
	//     }
	// 	for (int i=0;i<6;i++){
	// 		if (!isfinite(Fglobal[i])) {
	// 			cerr << "[ERROR] Fglobal NaN at elem " << this->elemID << " idx="<<i<<"\n";
	// 			Fglobal[i] = 0.0;
	// 		}
	// 	}
	//     // assemble into global vector fPtr->F
	//     assembleGlobalMat(Fglobal, fPtr->F);
	// }
// CORRECT
// void FEMElementClass::forceVector() {
//     // 获取当前几何
//     double xi = node[0]->pos[0], yi = node[0]->pos[1];
//     double xj = node[1]->pos[0], yj = node[1]->pos[1];
    
//     // 当前长度
//     L = sqrt((xj - xi)*(xj - xi) + (yj - yi)*(yj - yi));
    
//     // === 关键：共旋框架下的局部变形 ===
    
//     // 计算刚体旋转
//     double rigid_rotation = atan2(yj - yi, xj - xi);
    
//     // 计算局部坐标系下的纯变形
//     double local_axial_strain = (L - L0) / L0;  // 轴向应变
    
//     // 节点在局部坐标系下的转角（相对于刚体旋转）
//     double local_theta1 = node[0]->angle - rigid_rotation;
//     double local_theta2 = node[1]->angle - rigid_rotation;
    
//     // 确保角度在合理范围内
//     local_theta1 = Utils::shiftAngle(local_theta1);
//     local_theta2 = Utils::shiftAngle(local_theta2);
    
//     // === 在局部坐标系中计算内力 ===
    
//     // 轴向力（基于小应变，因为大变形已被刚体旋转吸收）
//     double F0 = (E * A / L0) * local_axial_strain * L0;  // 乘以L0得到力
    
//     // 弯矩（在局部坐标系中是线性的）
//     double M1 = (2 * E * I / L0) * (2.0 * local_theta1 + local_theta2);
//     double M2 = (2 * E * I / L0) * (local_theta1 + 2.0 * local_theta2);
    
//     // === 组装局部内力向量 ===
//     F[0] = -F0;
//     F[1] = (M1 + M2) / L0;  // 使用初始长度
//     F[2] = M1;
//     F[3] = F0;
//     F[4] = -(M1 + M2) / L0;
//     F[5] = M2;
    
//     // axialForce = F0;
//     // shearForce = (M1 + M2) / L0;
    
//     // 诊断
//     if (elemID == 0) {
//         cout << "Element " << elemID << " COROTATIONAL diagnostics:" << endl;
//         cout << "  Rigid rotation = " << rigid_rotation * 180.0 / M_PI << " deg" << endl;
//         cout << "  Local thetas = " << local_theta1 * 180.0 / M_PI << ", " 
//              << local_theta2 * 180.0 / M_PI << " deg" << endl;
//         cout << "  Axial strain = " << local_axial_strain << endl;
//     }
    
//     // 转换到全局坐标系
//     array<double, elementDOFs> FGlobal = Utils::Transpose(T) * F;
//     assembleGlobalMat(FGlobal, fPtr->F);
// }

// void FEMElementClass::forceVector() {
//     // 1. 当前几何
//     double xi = node[0]->pos[0];
//     double yi = node[0]->pos[1];
//     double xj = node[1]->pos[0];
//     double yj = node[1]->pos[1];

//     L = sqrt((xj - xi)*(xj - xi) + (yj - yi)*(yj - yi));
//     if (!isfinite(L) || L < 1e-12) L = L0;
    
//     c = (xj - xi) / L;
//     s = (yj - yi) / L;
//     if (!isfinite(c) || !isfinite(s)) { c = 1.0; s = 0.0; }

//     // 2. 局部变形
//     double theta1 = Utils::shiftAngle(node[0]->angle - angle);
//     double theta2 = Utils::shiftAngle(node[1]->angle - angle);

//     // === 3. 构造局部刚度矩阵（与 stiffMatrix 完全一致）===
//     double EA = E * A;
//     double EI = E * I;
//     double L02 = L0 * L0;
//     double L03 = L02 * L0;

//     double K_loc[6][6] = {0};
    
//     K_loc[0][0] =  EA / L0;   K_loc[0][3] = -EA / L0;
//     K_loc[3][0] = -EA / L0;   K_loc[3][3] =  EA / L0;
    
//     K_loc[1][1] = 12.0*EI/L03;   K_loc[1][2] = 6.0*EI/L02;
//     K_loc[1][4] = -12.0*EI/L03;  K_loc[1][5] = 6.0*EI/L02;
//     K_loc[2][1] = 6.0*EI/L02;    K_loc[2][2] = 4.0*EI/L0;
//     K_loc[2][4] = -6.0*EI/L02;   K_loc[2][5] = 2.0*EI/L0;
//     K_loc[4][1] = -12.0*EI/L03;  K_loc[4][2] = -6.0*EI/L02;
//     K_loc[4][4] = 12.0*EI/L03;   K_loc[4][5] = -6.0*EI/L02;
//     K_loc[5][1] = 6.0*EI/L02;    K_loc[5][2] = 2.0*EI/L0;
//     K_loc[5][4] = -6.0*EI/L02;   K_loc[5][5] = 4.0*EI/L0;

//     // === 4. 局部位移向量（关键修正！）===
//     // 在局部坐标系中，节点1在原点，节点2在 (L0, 0)
//     // 位移 = 当前位置 - 初始位置
    
//     // 对于 corotational 方法：
//     // - 轴向：两个节点都没有轴向位移（因为已经通过 L 的变化考虑了）
//     // - 横向：两个节点都没有横向位移（corotational 假设）
//     // - 转角：theta1 和 theta2 是相对单元的转角
    
//     double u_loc[6] = {
//         0.0,      // 节点1 轴向位移（局部坐标）
//         0.0,      // 节点1 横向位移（局部坐标）
//         theta1,   // 节点1 转角
//         0.0,      // 节点2 轴向位移（局部坐标）
//         0.0,      // 节点2 横向位移（局部坐标）
//         theta2    // 节点2 转角
//     };
    
//     // === 5. 内力 = K * u（只计算弯曲内力）===
//     for (int i = 0; i < 6; i++) {
//         Flocal[i] = 0.0;
//         for (int j = 0; j < 6; j++) {
//             Flocal[i] += K_loc[i][j] * u_loc[j];
//         }
//     }
    
//     // === 6. 添加轴向内力（单独计算）===
//     // 轴向力由长度变化产生
//     double u_axial = L - L0;
//     double N = (EA / L0) * u_axial;
    
//     Flocal[0] = -N;  // 节点1的轴向内力
//     Flocal[3] = N;   // 节点2的轴向内力

//     // === 7. 转换到全局坐标 ===
//     double Tmat[6][6] = {
//         { c,  s, 0,  0,  0, 0 },
//         { -s, c, 0,  0,  0, 0 },
//         { 0,  0, 1,  0,  0, 0 },
//         { 0,  0, 0,  c,  s, 0 },
//         { 0,  0, 0, -s,  c, 0 },
//         { 0,  0, 0,  0,  0, 1 }
//     };

//     std::array<double,6> Fglobal;
//     for (int i=0;i<6;i++){
//         Fglobal[i] = 0.0;
//         for (int j=0;j<6;j++){
//             Fglobal[i] += Tmat[j][i] * Flocal[j];
//         }
//     }
    
//     for (int i=0;i<6;i++){
//         if (!isfinite(Fglobal[i])) {
//             cerr << "[ERROR] Fglobal NaN at elem " << elemID << " idx="<<i<<"\n";
//             Fglobal[i] = 0.0;
//         }
//     }

//     assembleGlobalMat(Fglobal, fPtr->F);
    
//     // === 诊断 ===
//     static int diag_force = 0;
//     if (diag_force < 5 && elemID == 0) {
//         cout << "[forceVector] elem=" << elemID << endl;
//         cout << "  L=" << L << ", L0=" << L0 << ", N=" << N << " N" << endl;
//         cout << "  theta1=" << theta1*180/M_PI << "°, theta2=" << theta2*180/M_PI << "°" << endl;
//         cout << "  Flocal = [" << Flocal[0] << ", " << Flocal[1] << ", " 
//              << Flocal[2] << ", " << Flocal[3] << ", " 
//              << Flocal[4] << ", " << Flocal[5] << "]" << endl;
//         diag_force++;
//     }
// }

// Construct mass matrix
void FEMElementClass::massMatrix() {

	// Get linear stiffness matrix in global coordinates
	array<array<double, elementDOFs>, elementDOFs> MGlobal = Utils::Transpose(T) * M * T;

	// Assemble into global matrix
	assembleGlobalMat(MGlobal, fPtr->M);
}


// Construct stiffness matrix
void FEMElementClass::stiffMatrix() {

	// Get linear stiffness matrix in global coordinates
	array<array<double, elementDOFs>, elementDOFs> KGlobal = Utils::Transpose(T) * K_L * T;

	// Internal forces
	double F0 = -F[0];
	double V0 = F[4];

	// Construct upper half of local stiffness matrix for single element
	K_NL[0][1] = -V0 / L0;
	K_NL[0][4] = V0 / L0;
	K_NL[1][0] = -V0 / L0;
	K_NL[1][1] = F0 / L0;
	K_NL[1][3] = V0 / L0;
	K_NL[1][4] = -F0 / L0;
	K_NL[3][1] = V0 / L0;
	K_NL[3][4] = -V0 / L0;
	K_NL[4][0] = V0 / L0;
	K_NL[4][1] = -F0 / L0;
	K_NL[4][3] = -V0 / L0;
	K_NL[4][4] = F0 / L0;

	// Multiply by transformation matrices to get global matrix for single element
	KGlobal = KGlobal + Utils::Transpose(T) * K_NL * T;

	// Assemble into global matrix
	assembleGlobalMat(KGlobal, fPtr->K);
}

// void FEMElementClass::stiffMatrix() {
//     // 1. 更新当前几何
//     double xi = node[0]->pos[0];
//     double yi = node[0]->pos[1];
//     double xj = node[1]->pos[0];
//     double yj = node[1]->pos[1];

//     L = sqrt((xj - xi)*(xj - xi) + (yj - yi)*(yj - yi));
//     if (!isfinite(L) || L < 1e-12) {
//         L = L0;
//     }
//     c = (xj - xi) / L;
//     s = (yj - yi) / L;
//     if (!isfinite(c) || !isfinite(s)) {
//         c = 1.0; s = 0.0;
//     }

//     // 2. 局部线性刚度矩阵（基于 L0）
//     double EA = E * A;
//     double EI = E * I;
//     double L0_2 = L0 * L0;
//     double L0_3 = L0_2 * L0;

//     for (int i = 0; i < 6; i++)
//         for (int j = 0; j < 6; j++)
//             Klocal[i][j] = 0.0;

//     // 轴向
//     Klocal[0][0] =  EA / L0;
//     Klocal[0][3] = -EA / L0;
//     Klocal[3][0] = -EA / L0;
//     Klocal[3][3] =  EA / L0;

//     // 弯曲（基于 L0）
//     Klocal[1][1] = 12.0 * EI / L0_3;
//     Klocal[1][2] = 6.0 * EI / L0_2;
//     Klocal[1][4] = -12.0 * EI / L0_3;
//     Klocal[1][5] = 6.0 * EI / L0_2;
    
//     Klocal[2][1] = 6.0 * EI / L0_2;
//     Klocal[2][2] = 4.0 * EI / L0;
//     Klocal[2][4] = -6.0 * EI / L0_2;
//     Klocal[2][5] = 2.0 * EI / L0;
    
//     Klocal[4][1] = -12.0 * EI / L0_3;
//     Klocal[4][2] = -6.0 * EI / L0_2;
//     Klocal[4][4] = 12.0 * EI / L0_3;
//     Klocal[4][5] = -6.0 * EI / L0_2;
    
//     Klocal[5][1] = 6.0 * EI / L0_2;
//     Klocal[5][2] = 2.0 * EI / L0;
//     Klocal[5][4] = -6.0 * EI / L0_2;
//     Klocal[5][5] = 4.0 * EI / L0;

//     // *** 3. 几何刚度矩阵（关键！）***
//     // 从内力向量获取轴向力
//     double N = -F[0];  // 轴向力（拉为正）
    
//     std::array<std::array<double,6>,6> Kg;
//     for (int i = 0; i < 6; i++)
//         for (int j = 0; j < 6; j++)
//             Kg[i][j] = 0.0;
    
//     // 几何刚度（轴向力导致的横向刚度变化）
//     double coef = N / L0;
    
//     Kg[1][1] = coef * 6.0/5.0;
//     Kg[1][2] = coef * L0/10.0;
//     Kg[1][4] = -coef * 6.0/5.0;
//     Kg[1][5] = coef * L0/10.0;
    
//     Kg[2][1] = coef * L0/10.0;
//     Kg[2][2] = coef * 2.0*L0_2/15.0;
//     Kg[2][4] = -coef * L0/10.0;
//     Kg[2][5] = -coef * L0_2/30.0;
    
//     Kg[4][1] = -coef * 6.0/5.0;
//     Kg[4][2] = -coef * L0/10.0;
//     Kg[4][4] = coef * 6.0/5.0;
//     Kg[4][5] = -coef * L0/10.0;
    
//     Kg[5][1] = coef * L0/10.0;
//     Kg[5][2] = -coef * L0_2/30.0;
//     Kg[5][4] = -coef * L0/10.0;
//     Kg[5][5] = coef * 2.0*L0_2/15.0;

//     // 总局部刚度 = 线性刚度 + 几何刚度
//     for (int i = 0; i < 6; i++)
//         for (int j = 0; j < 6; j++)
//             Klocal[i][j] += Kg[i][j];

//     // 4. 旋转到全局坐标系
//     double Tmat[6][6] = {
//         { c,  s, 0,  0,  0, 0 },
//         {-s,  c, 0,  0,  0, 0 },
//         { 0,  0, 1,  0,  0, 0 },
//         { 0,  0, 0,  c,  s, 0 },
//         { 0,  0, 0, -s,  c, 0 },
//         { 0,  0, 0,  0,  0, 1 }
//     };

//     std::array<std::array<double,6>,6> temp;
//     for (int i = 0; i < 6; i++) {
//         for (int j = 0; j < 6; j++) {
//             temp[i][j] = 0.0;
//             for (int k = 0; k < 6; k++)
//                 temp[i][j] += Klocal[i][k] * Tmat[k][j];
//         }
//     }
    
//     for (int i = 0; i < 6; i++) {
//         for (int j = 0; j < 6; j++) {
//             Kglobal[i][j] = 0.0;
//             for (int k = 0; k < 6; k++)
//                 Kglobal[i][j] += Tmat[k][i] * temp[k][j];
//         }
//     }

//     // 5. 检查并组装
//     for (int i = 0; i < 6; i++) {
//         for (int j = 0; j < 6; j++) {
//             if (!isfinite(Kglobal[i][j])) {
//                 cerr << "[ERROR] Kglobal NaN at elem " << elemID << "\n";
//                 Kglobal[i][j] = 0.0;
//             }
//         }
//     }
// 	// 在stiffMatrix()中添加诊断
// 	if (elemID == 0 && fabs(N) > 1e-6) {
// 		cout << "Element " << elemID << ": N=" << N 
// 			<< ", Geometric stiffness factor=" << (N/L0) << endl;
// 	}
//     assembleGlobalMat(Kglobal, fPtr->K);
// }


//CORRECT ANGLE
// void FEMElementClass::stiffMatrix() {
//     // 1) current geometry (for transform) - used only to compute c,s to rotate local->global
//     double xi = node[0]->pos[0];
//     double yi = node[0]->pos[1];
//     double xj = node[1]->pos[0];
//     double yj = node[1]->pos[1];

//     L = sqrt((xj - xi)*(xj - xi) + (yj - yi)*(yj - yi));
//     if (!isfinite(L) || L < 1e-12) L = L0;

//     c = (xj - xi) / L;
//     s = (yj - yi) / L;
//     if (!isfinite(c) || !isfinite(s)) { c=1.0; s=0.0; }

//     // 2) linear local stiffness based on reference length L0 (consistent with reference force)
//     double EA = E * A;
//     double EI = E * I;
//     double L03 = L0 * L0 * L0;
//     double L02 = L0 * L0;

//     double Kloc[6][6] = {0};

//     // axial
//     Kloc[0][0] =  EA / L0;   Kloc[0][3] = -EA / L0;
//     Kloc[3][0] = -EA / L0;   Kloc[3][3] =  EA / L0;

//     // bending (use L0 here to be consistent with reference-based R)
//     Kloc[1][1] = 12.0*EI / L03;   Kloc[1][2] = 6.0*EI / L02;
//     Kloc[1][4] = -12.0*EI / L03;  Kloc[1][5] = 6.0*EI / L02;

//     Kloc[2][1] = 6.0*EI / L02;    Kloc[2][2] = 4.0*EI / L0;
//     Kloc[2][4] = -6.0*EI / L02;   Kloc[2][5] = 2.0*EI / L0;

//     Kloc[4][1] = -12.0*EI / L03;  Kloc[4][2] = -6.0*EI / L02;
//     Kloc[4][4] = 12.0*EI / L03;   Kloc[4][5] = -6.0*EI / L02;

//     Kloc[5][1] = 6.0*EI / L02;    Kloc[5][2] = 2.0*EI / L0;
//     Kloc[5][4] = -6.0*EI / L02;   Kloc[5][5] = 4.0*EI / L0;

//     // 3) geometric stiffness: derive from Flocal axial force.
//     // Ensure Flocal is up-to-date (forceVector already computed it).
//     // Convention: Flocal[0] should equal -F0 (internal axial force on left node).
//     // We want N positive for tension; if Flocal[0] stores -F0, then N = -Flocal[0].
//     double N = 0.0;
//     // prefer to use element's Flocal (computed in forceVector)
//     if (isfinite(Flocal[0])) {
//         N = - Flocal[0]; // common convention: Flocal[0] = -F0, so -Flocal[0] gives axial internal force (tension positive)
//     } else {
//         N = 0.0;
//     }

// 	// // 计算当前轴向应变
//     // double u = L - L0;
    
//     // // 计算轴向力
//     // double N = (E * A / L0) * u;
    
//     // // 也可以从节点位移计算（更准确）
//     // double theta1 = Utils::shiftAngle(node[0]->angle - angle);
//     // double theta2 = Utils::shiftAngle(node[1]->angle - angle);
    
//     // // 考虑弯曲引起的轴向缩短（可选，但推荐）
//     // double u_bending_correction = -L0 * (theta1*theta1 + theta2*theta2 + theta1*theta2) / 12.0;
//     // double u_total = u + u_bending_correction;
    
//     // N = (E * A / L0) * u_total;
    

//     // Kg local (from standard beam geometric stiffness terms)
//     double coef = N / L0;
//     double Kg[6][6] = {0};

//     // populate Kg (symmetric) - coefficients from standard formulation (consistent scale)
//     Kg[1][1] = coef * 6.0/5.0;
//     Kg[1][2] = coef * L0/10.0;
//     Kg[1][4] = -coef * 6.0/5.0;
//     Kg[1][5] = coef * L0/10.0;

//     Kg[2][1] = coef * L0/10.0;
//     Kg[2][2] = coef * 2.0*L0*L0/15.0;
//     Kg[2][4] = -coef * L0/10.0;
//     Kg[2][5] = -coef * L0*L0/30.0;

//     Kg[4][1] = -coef * 6.0/5.0;
//     Kg[4][2] = -coef * L0/10.0;
//     Kg[4][4] = coef * 6.0/5.0;
//     Kg[4][5] = -coef * L0/10.0;

//     Kg[5][1] = coef * L0/10.0;
//     Kg[5][2] = -coef * L0*L0/30.0;
//     Kg[5][4] = -coef * L0/10.0;
//     Kg[5][5] = coef * 2.0*L0*L0/15.0;

//     // Optionally scale Kg for diagnostics:
//     double alphaKg = 1.0; // try 0.0,0.5,1.0
//     for (int i=0;i<6;i++) for (int j=0;j<6;j++) Kloc[i][j] += alphaKg * Kg[i][j];

//     // 4) rotate to global: Kglob = T^T * Kloc * T where T built from current c,s
//     double T[6][6] = {
//         { c,  s, 0,  0,  0, 0 },
//         {-s,  c, 0,  0,  0, 0 },
//         { 0,   0, 1,  0,  0, 0 },
//         { 0,   0,  0,  c,  s, 0 },
//         { 0,   0,  0, -s,  c, 0 },
//         { 0,   0,  0,  0,  0, 1 }
//     };

//     double tmp[6][6] = {0}, Kglob[6][6] = {0};
//     for (int i=0;i<6;i++) for (int j=0;j<6;j++) {
//         for (int k=0;k<6;k++) tmp[i][j] += Kloc[i][k] * T[k][j];
//     }
//     for (int i=0;i<6;i++) for (int j=0;j<6;j++) {
//         for (int k=0;k<6;k++) Kglob[i][j] += T[k][i] * tmp[k][j];
//     }

//     // 5) assemble into global K with assembleGlobalMat
//     std::array<std::array<double, elementDOFs>, elementDOFs> Karr;
//     for (int i=0;i<6;i++) for (int j=0;j<6;j++) Karr[i][j] = Kglob[i][j];

//     // diagnostics
//     if (elemID == 0) {
//         cerr << "Element " << elemID << " diag: N="<<N<<" coef="<<coef
//              <<" alphaKg="<<alphaKg<<"\n";
//     }

//     // static int diag_count = 0;
//     // if (diag_count < 10 && elemID == 0) {
//     //     cout << "Element " << elemID << " stiffMatrix:" << endl;
//     //     cout << "  L=" << L << ", L0=" << L0 << ", u=" << u << endl;
//     //     cout << "  theta1=" << theta1*180/M_PI << "°, theta2=" << theta2*180/M_PI << "°" << endl;
//     //     cout << "  u_bending_correction=" << u_bending_correction << endl;
//     //     cout << "  N=" << N << " N (axial force)" << endl;
//     //     cout << "  coef=" << coef << endl;
//     //     diag_count++;
//     // }
//     assembleGlobalMat(Karr, fPtr->K);
// }

// void FEMElementClass::stiffMatrix() {
//     // 1. 当前几何
//     double xi = node[0]->pos[0];
//     double yi = node[0]->pos[1];
//     double xj = node[1]->pos[0];
//     double yj = node[1]->pos[1];

//     L = sqrt((xj - xi)*(xj - xi) + (yj - yi)*(yj - yi));
//     if (!isfinite(L) || L < 1e-12) L = L0;

//     c = (xj - xi) / L;
//     s = (yj - yi) / L;
//     if (!isfinite(c) || !isfinite(s)) { c=1.0; s=0.0; }

//     // 2. 线性局部刚度（基于 L0）
//     double EA = E * A;
//     double EI = E * I;
//     double L03 = L0 * L0 * L0;
//     double L02 = L0 * L0;

//     double Kloc[6][6] = {0};

//     Kloc[0][0] =  EA / L0;   Kloc[0][3] = -EA / L0;
//     Kloc[3][0] = -EA / L0;   Kloc[3][3] =  EA / L0;

//     Kloc[1][1] = 12.0*EI / L03;   Kloc[1][2] = 6.0*EI / L02;
//     Kloc[1][4] = -12.0*EI / L03;  Kloc[1][5] = 6.0*EI / L02;
//     Kloc[2][1] = 6.0*EI / L02;    Kloc[2][2] = 4.0*EI / L0;
//     Kloc[2][4] = -6.0*EI / L02;   Kloc[2][5] = 2.0*EI / L0;
//     Kloc[4][1] = -12.0*EI / L03;  Kloc[4][2] = -6.0*EI / L02;
//     Kloc[4][4] = 12.0*EI / L03;   Kloc[4][5] = -6.0*EI / L02;
//     Kloc[5][1] = 6.0*EI / L02;    Kloc[5][2] = 2.0*EI / L0;
//     Kloc[5][4] = -6.0*EI / L02;   Kloc[5][5] = 4.0*EI / L0;

//     // === 3. 几何刚度（关键修正）===
//     // 从当前变形计算轴向力
//     double u_axial = L - L0;
//     double N = (EA / L0) * u_axial;  // 正确的轴向力
    
//     // 限幅（防止数值问题）
//     double N_max = EA * 0.1;  // 10% 应变上限
//     if (N > N_max) N = N_max;
//     if (N < -N_max) N = -N_max;
    
//     double coef = N / L0;
//     double Kg[6][6] = {0};

//     Kg[1][1] = coef * 6.0/5.0;
//     Kg[1][2] = coef * L0/10.0;
//     Kg[1][4] = -coef * 6.0/5.0;
//     Kg[1][5] = coef * L0/10.0;

//     Kg[2][1] = coef * L0/10.0;
//     Kg[2][2] = coef * 2.0*L0*L0/15.0;
//     Kg[2][4] = -coef * L0/10.0;
//     Kg[2][5] = -coef * L0*L0/30.0;

//     Kg[4][1] = -coef * 6.0/5.0;
//     Kg[4][2] = -coef * L0/10.0;
//     Kg[4][4] = coef * 6.0/5.0;
//     Kg[4][5] = -coef * L0/10.0;

//     Kg[5][1] = coef * L0/10.0;
//     Kg[5][2] = -coef * L0*L0/30.0;
//     Kg[5][4] = -coef * L0/10.0;
//     Kg[5][5] = coef * 2.0*L0*L0/15.0;

//     // 总刚度
//     for (int i=0;i<6;i++) 
//         for (int j=0;j<6;j++) 
//             Kloc[i][j] += Kg[i][j];

//     // 4. 转换到全局
//     double T[6][6] = {
//         { c,  s, 0,  0,  0, 0 },
//         {-s,  c, 0,  0,  0, 0 },
//         { 0,  0, 1,  0,  0, 0 },
//         { 0,  0, 0,  c,  s, 0 },
//         { 0,  0, 0, -s,  c, 0 },
//         { 0,  0, 0,  0,  0, 1 }
//     };

//     double tmp[6][6] = {0}, Kglob[6][6] = {0};
//     for (int i=0;i<6;i++) for (int j=0;j<6;j++) {
//         for (int k=0;k<6;k++) tmp[i][j] += Kloc[i][k] * T[k][j];
//     }
//     for (int i=0;i<6;i++) for (int j=0;j<6;j++) {
//         for (int k=0;k<6;k++) Kglob[i][j] += T[k][i] * tmp[k][j];
//     }

//     // 5. 组装
//     std::array<std::array<double, elementDOFs>, elementDOFs> Karr;
//     for (int i=0;i<6;i++) for (int j=0;j<6;j++) Karr[i][j] = Kglob[i][j];

//     // 诊断
//     static int diag_count = 0;
//     if (diag_count < 5 && elemID == 0) {
//         cout << "Element " << elemID << " stiffMatrix:" << endl;
//         cout << "  L=" << L << ", L0=" << L0 << ", u=" << u_axial << endl;
//         cout << "  N=" << N << " N, coef=" << coef << endl;
//         diag_count++;
//     }

//     assembleGlobalMat(Karr, fPtr->K);
// }

// Multiply by shape functions
array<double, dims> FEMElementClass::shapeFuns(const array<double, elementDOFs> &vec, double zeta) {

	// Results vector
	array<double, dims> resVec;

	// Use shape functions to calculate values
	double N0 = 1.0 - (zeta + 1.0) / 2.0;
	double N1 = 1.0 - 3.0 * SQ((zeta + 1.0) / 2.0) + 2.0 * TH((zeta + 1.0) / 2.0);
	double N2 = ((zeta + 1.0) / 2.0 - 2.0 * SQ((zeta + 1.0) / 2.0) + TH((zeta + 1.0) / 2.0)) * L;
	double N3 = (zeta + 1.0) / 2.0;
	double N4 = 3.0 * SQ((zeta + 1.0) / 2.0) - 2.0 * TH((zeta + 1.0) / 2.0);
	double N5 = (-SQ((zeta + 1.0) / 2.0) + TH((zeta + 1.0) / 2.0)) * L;

	// Calculate values using shape functions
	resVec[eX] = vec[0] * N0 + vec[3] * N3;
	resVec[eY] = vec[1] * N1 + vec[2] * N2 + vec[4] * N4 + vec[5] * N5;

	// Return
	return resVec;
}

// Assemble into global vector
void FEMElementClass::assembleGlobalMat(const array<double, elementDOFs> &localVec, vector<double> &globalVec) {

	// Loop through and set
	for (size_t i = 0; i < localVec.size(); i++)
		globalVec[DOFs[i]] += localVec[i];
}

// Assemble into global matrix
void FEMElementClass::assembleGlobalMat(const array<array<double, elementDOFs>, elementDOFs> &localVec, vector<double> &globalVec) {

	// Get dimensions
	int dim = fPtr->bodyDOFs;

	// Get rows and cols
	size_t rows = localVec.size();
	size_t cols = localVec[0].size();

	// Now loop through and set
	for (size_t i = 0; i < rows; i++) {
		for (size_t j = 0; j < cols; j++) {
			globalVec[DOFs[i] * dim + DOFs[j]] += localVec[i][j];
		}
	}
}

// Assemble into global vector
array<double, elementDOFs> FEMElementClass::disassembleGlobalMat(const vector<double> &globalVec) {

	// Declare vector
	array<double, elementDOFs> localVec;

	// Loop through and set
	for (size_t i = 0; i < localVec.size(); i++)
		localVec[i] = globalVec[DOFs[i]];

	// Return
	return localVec;
}

// Set local mass and stiffness matrices
void FEMElementClass::setLocalMatrices() {

	// Size the local matrices
	M.fill({0.0});
	K_L.fill({0.0});
	K_NL.fill({0.0});
	R.fill(0.0);
	F.fill(0.0);

	// Coefficients
	double C1 = rho * A * L0 / 420.0;

	// Set mass matrix first
	M[0][0] = C1 * 140.0;
	M[0][3] = C1 * 70.0;
	M[1][1] = C1 * 156.0;
	M[1][2] = C1 * 22.0 * L0;
	M[1][4] = C1 * 54;
	M[1][5] = C1 * (-13.0 * L0);
	M[2][2] = C1 * 4.0 * SQ(L0);
	M[2][4] = C1 * 13.0 * L0;
	M[2][5] = C1 * (-3.0 * SQ(L0));
	M[3][3] = C1 * 140.0;
	M[4][4] = C1 * 156.0;
	M[4][5] = C1 * (-22.0 * L0);
	M[5][5] = C1 * 4.0 * SQ(L0);

	// Now set stiffness matrix
	K_L[0][0] = E * A / L0;
	K_L[0][3] = -E * A / L0;
	K_L[1][1] = 12.0 * E * I / TH(L0);
	K_L[1][2] = 6.0 * E * I / SQ(L0);
	K_L[1][4] = -12.0 * E * I / TH(L0);
	K_L[1][5] = 6.0 * E * I / SQ(L0);
	K_L[2][2] = 4.0 * E * I / L0;
	K_L[2][4] = -6.0 * E * I / SQ(L0);
	K_L[2][5] = 2.0 * E * I / L0;
	K_L[3][3] = E * A / L0;
	K_L[4][4] = 12.0 * E * I / TH(L0);
	K_L[4][5] = -6.0 * E * I / SQ(L0);
	K_L[5][5] = 4.0 * E * I / L0;

	// Copy to the lower half (symmetrical matrix)
	for (size_t i = 1; i < M.size(); i++) {
		for (size_t j = 0; j < i; j++) {
			M[i][j] = M[j][i];
			K_L[i][j] = K_L[j][i];
		}
	}
}

// Set transformation matrix for element
void FEMElementClass::setElementTransform() {

	// Set to correct values
	T[0][0] = T[1][1] =  T[3][3] = T[4][4] = cos(angle);
	T[0][1] = T[3][4] = sin(angle);
	T[1][0] = T[4][3] = -sin(angle);
	T[2][2] = T[5][5] =  1.0;
}

// Custom constructor for building elements
FEMElementClass::FEMElementClass(FEMBodyClass *fBody, int i, const array<double, dims> &geom, double angleRad, double length, double den, double youngMod) {

	// Set pointer to fBody
	fPtr = fBody;

	// Set values
	L0 = length;
	//L = length;
	angle = angleRad;
	// double b = fPtr->iPtr->oPtr->gPtr->Dx;   // physical length per lattice spacing
    // double h = geom[1];                      // height (as provided in geom) -- make sure this is in same units (physical or LU)
	// A = b * h;
    // I = b * h * h * h / 12.0;  
	A = fPtr->iPtr->oPtr->gPtr->Dx * geom[1];
	I = fPtr->iPtr->oPtr->gPtr->Dx * TH(geom[1]) / 12.0;
	E = youngMod;
	rho = den;
	// double D = geom[1];  
    
    // // 圆形截面（论文用的）
    // A = M_PI * D * D / 4.0;           // πD²/4
    // I = M_PI * D * D * D * D / 64.0;  // πD⁴/64
    

	// Vector of node indices
	array<int, elementNodes> nodeIdx;

	// Get nodal indices
	nodeIdx[0] = i;
	nodeIdx[1] = i + 1;

	// Set pointers and global DOFs
	for (int i = 0; i < elementNodes; i++) {

		// Insert pointer to node
		node[i] = &(fPtr->node[nodeIdx[i]]);

		// Get DOFs
		for (int j = 0; j < nodeDOFs; j++)
			DOFs[i * nodeDOFs + j] = nodeIdx[i] * nodeDOFs + j;
	}

	// Set transformation matrix
	T.fill({0.0});
	setElementTransform();

	// Set local mass and stiffness matrices
	setLocalMatrices();
}
