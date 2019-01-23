	Eigen::MatrixXd Vout = V;
		Eigen::MatrixXi Fout = F;

		// Add new node at triangle center
		igl::adjacency_list(F, VV, true);
		float divide = 0.3334;
		Vout.conservativeResize((V.rows() + F.rows()), 3);
		Fout.conservativeResize((F.rows() + F.rows() *2), 3);

		for (int f = 0; f < F.rows(); ++f) {
			auto v1i = F(f, 0);    auto v1 = V.row(v1i);
			auto v2i = F(f, 1);    auto v2 = V.row(v2i);
			auto v3i = F(f, 2);    auto v3 = V.row(v3i);
			
			int insertVIndex = V.rows()+ f;
			Vout.row(insertVIndex)  <<	divide * (v1(0) + v2(0) + v3(0)),
									divide * (v1(1) + v2(1) + v3(1)),
									divide * (v1(2) + v2(2) + v3(2));
			
			Fout.row(f)<<					v1i, v2i, insertVIndex;
			Fout.row(F.rows() + 2* f) <<		v2i, v3i, insertVIndex;
			Fout.row(F.rows() + 2* f+1) <<	v3i, v1i, insertVIndex;

		}

		// Move origin vertices
		for (int i = 0; i < V.rows(); ++i) {
			int n = VV[i].size();
			auto a_n = (4 - 2 * cos(M_PI * 2 / n)) / 9;
			Eigen::Vector3d  sum = { 0, 0, 0 };

			for (auto nb : VV[i]) {
				sum += V.row(nb);
			}
		
			Eigen::Vector3d p =V.row(i);
			auto newP = p* (1 - a_n) + a_n / n * sum;
			Vout.row(i) << newP(0), newP(1), newP(2);
		}
		
		// Flip edge
		igl::vertex_triangle_adjacency(Vout, Fout, VF, VFi);
		for (int ov1 = 0; ov1 < V.rows(); ov1++) {
			auto nb_ov1 = VV[ov1];
			for (auto ov2 : nb_ov1) {
				if (ov2 < ov1) {
					continue;
				}
				// flip  edge of ov1 to ov2
				std::vector<int> nb_faces;
				for (auto f : VF[ov1]) {
					if (std::find(VF[ov2].begin(), VF[ov2].end(), f) != VF[ov2].end()) {
						nb_faces.push_back(f);
					}
				}
				//flip two face
				if (nb_faces.size() != 2) {
					cout << "error" << endl;		continue;
				}

				std::vector<int> new_verts;
				for (auto fi : nb_faces) {
					auto f = Fout.row(fi);
					for (int i = 0; i < 3; i++) {
						if (f(i) != ov1&& f(i) != ov2) {
							new_verts.push_back(f(i));
						}
					}
				}
				
				if (new_verts.size() != 2) {
					cout << "error" << endl;    continue;
				}
				
				int fi, vi;
				Eigen::Vector3i f;
				int usedOldNode =-1;
				
				for (int i = 0; i < 2; i++) {
					fi= nb_faces[i];
					f = Fout.row(fi);
					vi = new_verts[1-i];
					
					for (int j = 0; j < 2; j++) {
						if (f(j) >= V.rows() && f(j + 1) < V.rows() && f(j + 1) != usedOldNode) {
							Fout.row(fi) << f(j), f(j + 1), vi;
							usedOldNode = f(j + 1);
							break;
						}
						else if (f(j + 1) >= V.rows() && f(j) < V.rows() && f(j) != usedOldNode) {
							Fout.row(fi) << vi, f(j), f(j + 1);
							usedOldNode = f(j);
							break;
						}
					}

				}
			}
		}

		V = Vout; F = Fout;
		update_display(viewer);
