#include <stdio.h>
#include <unordered_map>
#include <vector>
#include <complex>
#include <algorithm>

using complex = std::complex<double>;

#include <queue>
std::vector<double> cons;

bool close2(double a, double b) {
	return std::abs(a - b) < 1.0e-7;
}

int constant(double value) {
	for (int i = 0; i < cons.size(); i++) {
		if (close2(value, cons[i])) {
			return i;
		}
	}
	cons.push_back(value);
	return cons.size() - 1;
}

int exp_sz(int P) {
	return (P + 1) * (P + 1);
}

void flush_constants(FILE* fp) {
	fprintf(fp, "\n");
	fprintf(fp, "               .align         32\n");
	for (int i = 0; i < cons.size(); i++) {
		fprintf(fp, "C%03i:          .float         %-18.10e\n", i, cons[i]);
		fprintf(fp, "               .float         %-18.10e\n", cons[i]);
		fprintf(fp, "               .float         %-18.10e\n", cons[i]);
		fprintf(fp, "               .float         %-18.10e\n", cons[i]);
		fprintf(fp, "               .float         %-18.10e\n", cons[i]);
		fprintf(fp, "               .float         %-18.10e\n", cons[i]);
		fprintf(fp, "               .float         %-18.10e\n", cons[i]);
		fprintf(fp, "               .float         %-18.10e\n", cons[i]);
	}
	cons.clear();
}

int lindex(int l, int m) {
	return l * (l + 1) + m;
}

constexpr int cindex(int l, int m) {
	return l * (l + 1) / 2 + m;
}

std::vector<complex> spherical_singular_harmonic(int P, double x, double y, double z) {
	const double r2 = x * x + y * y + z * z;
	const double r2inv = double(1) / r2;
	complex R = complex(x, y);
	std::vector<complex> O(exp_sz(P));
	O[cindex(0, 0)] = complex(sqrt(r2inv), double(0));
	R *= r2inv;
	z *= r2inv;
	for (int m = 0; m <= P; m++) {
		if (m > 0) {
			O[cindex(m, m)] = O[cindex(m - 1, m - 1)] * R * double(2 * m - 1);
		}
		if (m + 1 <= P) {
			O[cindex(m + 1, m)] = double(2 * m + 1) * z * O[cindex(m, m)];
		}
		for (int n = m + 2; n <= P; n++) {
			O[cindex(n, m)] = (double(2 * n - 1) * z * O[cindex(n - 1, m)] - double((n - 1) * (n - 1) - m * m) * r2inv * O[cindex(n - 2, m)]);
		}
	}
	return O;
}

struct hash {
	size_t operator()(std::array<int, 3> i) const {
		return i[0] * 12345 + i[1] * 42 + i[2];
	}
};

double brot(int n, int m, int l) {
	static std::unordered_map<std::array<int, 3>, double, hash> values;
	std::array<int, 3> key;
	key[0] = n;
	key[1] = m;
	key[2] = l;
	if (values.find(key) != values.end()) {
		return values[key];
	} else {
		double v;
		if (n == 0 && m == 0 && l == 0) {
			v = 1.0;
		} else if (abs(l) > n) {
			v = 0.0;
		} else if (m == 0) {
			v = 0.5 * (brot(n - 1, m, l - 1) - brot(n - 1, m, l + 1));
		} else if (m > 0) {
			v = 0.5 * (brot(n - 1, m - 1, l - 1) + brot(n - 1, m - 1, l + 1) + 2.0 * brot(n - 1, m - 1, l));
		} else {
			v = 0.5 * (brot(n - 1, m + 1, l - 1) + brot(n - 1, m + 1, l + 1) - 2.0 * brot(n - 1, m + 1, l));
		}
		values[key] = v;
		return v;
	}
}

void greens(FILE* fp, int P) {
	fprintf(fp, "#define        x              %%ymm0\n");
	fprintf(fp, "#define        y              %%ymm1\n");
	fprintf(fp, "#define        z              %%ymm2\n");
	fprintf(fp, "#define        r2             %%ymm3\n");
	fprintf(fp, "#define        r2inv          %%ymm4\n");
	fprintf(fp, "#define        ax0            %%ymm5\n");
	fprintf(fp, "#define        ay0            %%ymm6\n");
	fprintf(fp, "#define        Optr           %%rdi\n");
	fprintf(fp, "\n");
	fprintf(fp, "               .global        greens\n");
	fprintf(fp, "\n");
	fprintf(fp, "               .text\n");
	fprintf(fp, "\n");
	fprintf(fp, "greens:        vmulps         z, z, r2\n");
	fprintf(fp, "               vfmadd231ps    y, y, r2\n");
	fprintf(fp, "               vfmadd231ps    x, x, r2\n");
	fprintf(fp, "               vmovaps        C%03i, %%ymm8\n", constant(-1.0));
	fprintf(fp, "               vdivps         r2, %%ymm8, r2inv\n");
	fprintf(fp, "               vrsqrtps       r2, %%ymm8\n");
	fprintf(fp, "               vmulps         C%03i, %%ymm8, %%ymm8\n", constant(-1.0));
	fprintf(fp, "               vmovaps        %%ymm8, (Optr)\n");
	fprintf(fp, "               vmulps         C%03i, r2inv, %%ymm9\n", constant(-1.0));
	fprintf(fp, "               vmulps         %%ymm9, x, x\n");
	fprintf(fp, "               vmulps         %%ymm9, y, y\n");
	fprintf(fp, "               vmulps         %%ymm9, z, z\n");
	fprintf(fp, "               vmulps         %%ymm8, x, %%ymm9\n");
	fprintf(fp, "               vmulps         %%ymm8, y, %%ymm8\n");
	fprintf(fp, "               vmovaps        %%ymm9, %i(Optr)\n", 32 * lindex(1, 1));
	fprintf(fp, "               vmovaps        %%ymm8, %i(Optr)\n", 32 * lindex(1, -1));
	for (int m = 2; m <= P; m++) {
		fprintf(fp, "               vmulps         C%03i, %%ymm9, ax0\n", constant(2 * m - 1));
		fprintf(fp, "               vmulps         C%03i, %%ymm8, ay0\n", constant(2 * m - 1));
		fprintf(fp, "               vmulps         y, ay0, %%ymm9\n");
		fprintf(fp, "               vmulps         x, ay0, %%ymm8\n");
		fprintf(fp, "               vfmsub231ps    x, ax0, %%ymm9\n");
		fprintf(fp, "               vfmadd231ps    y, ax0, %%ymm8\n");
		fprintf(fp, "               vmovaps        %%ymm9, %i(Optr)\n", 32 * lindex(m, m));
		fprintf(fp, "               vmovaps        %%ymm8, %i(Optr)\n", 32 * lindex(m, -m));
	}
	fprintf(fp, "               vmulps         (Optr), r2inv, %%ymm8\n");
	fprintf(fp, "               vmulps         C%03i, %%ymm8, %%ymm8\n", constant(-1.0));
	fprintf(fp, "               vmovaps        %%ymm8, %i(Optr)\n", 32 * lindex(2, 0));
	fprintf(fp, "               vmovaps        %%ymm9, %i(Optr)\n", 32 * lindex(1, 0));
	for (int n = 1; n < P; n++) {
		int mmax = n - 1;
		const double c0 = -(double((n + 1) * (n + 1)));
		fprintf(fp, "               vmulps         C%03i, z, ax0\n", constant(2 * n + 1));
		fprintf(fp, "               vmulps         %i(Optr), ax0, %%ymm8\n", 32 * lindex(n, -n));
		for (int m = -mmax; m <= mmax; m++) {
			fprintf(fp, "               vfmadd231ps    %i(Optr), ax0, %%ymm8\n", 32 * lindex(n, m));
		}
		fprintf(fp, "               vfmadd231ps    %i(Optr), ax0, %%ymm8\n", 32 * lindex(n, n));
		fprintf(fp, "               vmovaps        %%ymm8, %i(Optr)\n", 32 * lindex(n + 1, -n));
		if (n != P - 1) {
			fprintf(fp, "               vmulps         C%03i, r2inv, %%ymm8\n", constant(c0));
			fprintf(fp, "               vmulps         %i(Optr), %%ymm8, %%ymm8\n", 32 * lindex(n, 0));
			fprintf(fp, "               vmovaps        %%ymm8, %i(Optr)\n", 32 * lindex(n + 2, 0));
			for (int m = 1; m <= n; m++) {
				const double c0 = -(double((n + 1) * (n + 1)) - double(m * m));
				fprintf(fp, "               vmulps         C%03i, r2inv, ax0\n", constant(c0));
				fprintf(fp, "               vmulps         %i(Optr), ax0, %%ymm9\n", 32 * lindex(n, m));
				fprintf(fp, "               vmulps         %i(Optr), ax0, %%ymm8\n", 32 * lindex(n, -m));
				fprintf(fp, "               vmovaps        %%ymm9, %i(Optr)\n", 32 * lindex(n + 2, m));
				fprintf(fp, "               vmovaps        %%ymm8, %i(Optr)\n", 32 * lindex(n + 2, -m));
			}
		}
	}
	fprintf(fp, "               ret\n");
	flush_constants(fp);
}

void regular_harmonic(FILE* fp, int P) {
	fprintf(fp, "#define        x              %%ymm0\n");
	fprintf(fp, "#define        y              %%ymm1\n");
	fprintf(fp, "#define        z              %%ymm2\n");
	fprintf(fp, "#define        r2             %%ymm3\n");
	fprintf(fp, "#define        ax0            %%ymm4\n");
	fprintf(fp, "#define        ay0            %%ymm5\n");
	fprintf(fp, "#define        Y              %%rdi\n");
	fprintf(fp, "\n");
	fprintf(fp, "               .global        regular_harmonic\n");
	fprintf(fp, "\n");
	fprintf(fp, "               .text\n");
	fprintf(fp, "\n");
	fprintf(fp, "regular_harmonic:\n");
	//tprint("Y[0] = TCAST(1);\n");
	fprintf(fp, "               vmovaps        C%03i, %%ymm9\n", constant(1));
	fprintf(fp, "               vmovaps        %%ymm9, (Y)\n");
	fprintf(fp, "               vmovaps        %%ymm9, %%ymm8\n");
	for (int m = 1; m <= P; m++) {
		if (m - 1 > 0) {
			//tprint("ax0 = Y[%i] * TCAST(%.20e);\n", lindex(m - 1, m - 1), 1.0 / (2.0 * m));
			fprintf(fp, "               vmulps         C%03i, %%ymm8, ax0\n", constant(1.0 / (2.0 * m)));
			//tprint("ay0 = Y[%i] * TCAST(%.20e);\n", lindex(m - 1, -(m - 1)), 1.0 / (2.0 * m));
			fprintf(fp, "               vmulps         C%03i, %%ymm9, ay0\n", constant(1.0 / (2.0 * m)));
			//tprint("Y[%i] = x * ax0 - y * ay0;\n", lindex(m, m));
			fprintf(fp, "               vmulps         y, ay0, %%ymm8\n");
			fprintf(fp, "               vmulps         x, ay0, %%ymm9\n");
			fprintf(fp, "               vfmsub231ps    x, ax0, %%ymm8\n");
			fprintf(fp, "               vfmadd231ps    y, ax0, %%ymm9\n");
			fprintf(fp, "               vmovaps        %%ymm8, %i(Y)\n", 32 * lindex(m, m));
			fprintf(fp, "               vmovaps        %%ymm9, %i(Y)\n", 32 * lindex(m, -m));
			//tprint("Y[%i] = fma(y, ax0, x * ay0);\n", lindex(m, -m));
		} else {
			//	tprint("Y[%i] = x * TCAST(%.20e);\n", lindex(m, m), 1.0 / (2.0 * m));
			fprintf(fp, "               vmulps         C%03i, x, %%ymm8\n", constant(1.0 / (2.0 * m)));
			fprintf(fp, "               vmulps         C%03i, x, %%ymm9\n", constant(1.0 / (2.0 * m)));
			fprintf(fp, "               vmovaps        %%ymm8, %i(Y)\n", 32 * lindex(m, m));
			//	tprint("Y[%i] = y * TCAST(%.20e);\n", lindex(m, -m), 1.0 / (2.0 * m));
			fprintf(fp, "               vmovaps        %%ymm9, %i(Y)\n", 32 * lindex(m, -m));
		}
	}
	const double c0 = -0.25;
	const double c1 = double(1) / (double((1) * (1)));
	if (2 <= P) {
		//tprint("Y[%i] = TCAST(-0.25) * r2;\n", lindex(2, 0));
		fprintf(fp, "               vmulps         C%03i, r2, %%ymm9\n", constant(-0.25));
		fprintf(fp, "               vmovaps        %%ymm9, %i(Y)\n", 32 * lindex(2, 0));
	}
	if (1 <= P) {
		//tprint("Y[%i] = z;\n", lindex(1, 0));
		fprintf(fp, "               vmovaps        z, %i(Y)\n", 32 * lindex(1, 0));
	}
	for (int n = 1; n < P; n++) {
		const double c0 = -double(1) / (double((n + 2) * (n + 2)));
		const double c1 = double(2 * n + 1) / (double((n + 1) * (n + 1)));
		if (n + 2 <= P) {
			//	tprint("Y[%i] = TCAST(%.20e) * r2 * Y[%i];\n", lindex(n + 2, 0), c0, lindex(n, 0));
			fprintf(fp, "               vmulps         C%03i, r2, %%ymm9\n", constant(c0));
			fprintf(fp, "               vmulps         %i(Y), %%ymm9, %%ymm9\n", 32 * lindex(n, 0));
			fprintf(fp, "               vmovaps        %%ymm9, %i(Y)\n", 32 * lindex(n + 2, 0));
		}
		//tprint("Y[%i] = fma(TCAST(%.20e) * z, Y[%i], Y[%i]);\n", lindex(n + 1, 0), c1, lindex(n, 0), lindex(n + 1, 0));
		for (int m = 1; m <= n; m++) {
			const double c0 = -double(1) / (double((n + 2) * (n + 2)) - double(m * m));
			const double c1 = double(2 * n + 1) / (double((n + 1) * (n + 1)) - double(m * m));
			if (n + 2 <= P) {
				//tprint("ax0 = TCAST(%.20e) * r2;\n", c0);
				fprintf(fp, "               vmulps         C%03i, r2, ax0\n", constant(c0));
			}
			if (n + 2 <= P) {
				//	tprint("Y[%i] = ax0 * Y[%i];\n", lindex(n + 2, m), lindex(n, m));
				fprintf(fp, "               vmulps         %i(Y), ax0, %%ymm8\n", 32 * lindex(n, m));
				//	tprint("Y[%i] = ax0 * Y[%i];\n", lindex(n + 2, -m), lindex(n, -m));
				fprintf(fp, "               vmulps         %i(Y), ax0, %%ymm9\n", lindex(n, -m));
				fprintf(fp, "               vmovaps        %%ymm8, %i(Y)\n", 32 * lindex(n + 2, m));
				fprintf(fp, "               vmovaps        %%ymm9, %i(Y)\n", 32 * lindex(n + 2, -m));
			}
			if (n == m) {
				//	tprint("Y[%i] = z * Y[%i];\n", lindex(n + 1, m), lindex(n, m));
				fprintf(fp, "               vmulps         %i(Y), z, %%ymm8\n", 32 * lindex(n, m));
				//	tprint("Y[%i] = z * Y[%i];\n", lindex(n + 1, -m), lindex(n, -m));
				fprintf(fp, "               vmulps         %i(Y), z, %%ymm9\n", 32 * lindex(n, -m));
				fprintf(fp, "               vmovaps        %%ymm8, %i(Y)\n", 32 * lindex(n + 2, m));
				fprintf(fp, "               vmovaps        %%ymm9, %i(Y)\n", 32 * lindex(n + 2, -m));
			} else {
				//	tprint("ay0 = TCAST(%.20e) * z;\n", c1);
				fprintf(fp, "               vmulps         C%03i, z, ay0\n", constant(c1));
				//	tprint("Y[%i] = fma(ay0, Y[%i], Y[%i]);\n", lindex(n + 1, m), lindex(n, m), lindex(n + 1, m));
				fprintf(fp, "               vmovaps        %i(Y), %%ymm8\n", 32 * lindex(n + 1, m));
				fprintf(fp, "               vmovaps        %i(Y), %%ymm9\n", 32 * lindex(n + 1, -m));
				fprintf(fp, "               vfmadd231ps    %i(Y), ay0, %%ymm8\n", 32 * lindex(n, m));
				fprintf(fp, "               vfmadd231ps    %i(Y), ay0, %%ymm9\n", 32 * lindex(n, -m));
				fprintf(fp, "               vmovaps        %%ymm8, %i(Y)\n", 32 * lindex(n + 1, m));
				fprintf(fp, "               vmovaps        %%ymm9, %i(Y)\n", 32 * lindex(n + 1, -m));
				//	tprint("Y[%i] = fma(ay0, Y[%i], Y[%i]);\n", lindex(n + 1, -m), lindex(n, -m), lindex(n + 1, -m));
			}
		}
	}
	fprintf(fp, "               ret\n");
	flush_constants(fp);
}

void P2M(FILE* fp, int P) {
	fprintf(fp, "#define        x              %%ymm0\n");
	fprintf(fp, "#define        y              %%ymm1\n");
	fprintf(fp, "#define        z              %%ymm2\n");
	fprintf(fp, "#define        r2             %%ymm3\n");
	fprintf(fp, "#define        ax0            %%ymm4\n");
	fprintf(fp, "#define        ay0            %%ymm5\n");
	fprintf(fp, "#define        M              %%rdi\n");
	fprintf(fp, "#define        Y              %%rsi\n");
	fprintf(fp, "\n");
	fprintf(fp, "               .global        M2M\n");
	fprintf(fp, "\n");
	fprintf(fp, "               .text\n");
	fprintf(fp, "\n");
	fprintf(fp, "M2M:           push           %%rbp\n");
	fprintf(fp, "               mov            %%rsp, %%rbp\n");
	fprintf(fp, "               sub            $%i, %%rsp\n", 32 * exp_sz(P));
	fprintf(fp, "               and            $0xffffffffffffffc0, %%rsp\n");
	fprintf(fp, "               mov            %%rsp, %%rsi\n");
	//tprint("Y[0] = TCAST(1);\n");
	fprintf(fp, "               vmovaps        C%03i, %%ymm9\n", constant(1));
	fprintf(fp, "               vmovaps        %%ymm9, (Y)\n");
	fprintf(fp, "               vmovaps        %%ymm9, %%ymm8\n");
	for (int m = 1; m <= P; m++) {
		if (m - 1 > 0) {
			//tprint("ax0 = Y[%i] * TCAST(%.20e);\n", lindex(m - 1, m - 1), 1.0 / (2.0 * m));
			fprintf(fp, "               vmulps         C%03i, %%ymm8, ax0\n", constant(1.0 / (2.0 * m)));
			//tprint("ay0 = Y[%i] * TCAST(%.20e);\n", lindex(m - 1, -(m - 1)), 1.0 / (2.0 * m));
			fprintf(fp, "               vmulps         C%03i, %%ymm9, ay0\n", constant(1.0 / (2.0 * m)));
			//tprint("Y[%i] = x * ax0 - y * ay0;\n", lindex(m, m));
			fprintf(fp, "               vmulps         y, ay0, %%ymm8\n");
			fprintf(fp, "               vmulps         x, ay0, %%ymm9\n");
			fprintf(fp, "               vfmsub231ps    x, ax0, %%ymm8\n");
			fprintf(fp, "               vfmadd231ps    y, ax0, %%ymm9\n");
			fprintf(fp, "               vmovaps        %%ymm8, %i(Y)\n", 32 * lindex(m, m));
			fprintf(fp, "               vmovaps        %%ymm9, %i(Y)\n", 32 * lindex(m, -m));
			//tprint("Y[%i] = fma(y, ax0, x * ay0);\n", lindex(m, -m));
		} else {
			//	tprint("Y[%i] = x * TCAST(%.20e);\n", lindex(m, m), 1.0 / (2.0 * m));
			fprintf(fp, "               vmulps         C%03i, x, %%ymm8\n", constant(1.0 / (2.0 * m)));
			fprintf(fp, "               vmulps         C%03i, x, %%ymm9\n", constant(1.0 / (2.0 * m)));
			fprintf(fp, "               vmovaps        %%ymm8, %i(Y)\n", 32 * lindex(m, m));
			//	tprint("Y[%i] = y * TCAST(%.20e);\n", lindex(m, -m), 1.0 / (2.0 * m));
			fprintf(fp, "               vmovaps        %%ymm9, %i(Y)\n", 32 * lindex(m, -m));
		}
	}
	const double c0 = -0.25;
	const double c1 = double(1) / (double((1) * (1)));
	if (2 <= P) {
		//tprint("Y[%i] = TCAST(-0.25) * r2;\n", lindex(2, 0));
		fprintf(fp, "               vmulps         C%03i, r2, %%ymm9\n", constant(-0.25));
		fprintf(fp, "               vmovaps        %%ymm9, %i(Y)\n", 32 * lindex(2, 0));
	}
	if (1 <= P) {
		//tprint("Y[%i] = z;\n", lindex(1, 0));
		fprintf(fp, "               vmovaps        z, %i(Y)\n", 32 * lindex(1, 0));
	}
	for (int n = 1; n < P; n++) {
		const double c0 = -double(1) / (double((n + 2) * (n + 2)));
		const double c1 = double(2 * n + 1) / (double((n + 1) * (n + 1)));
		if (n + 2 <= P) {
			//	tprint("Y[%i] = TCAST(%.20e) * r2 * Y[%i];\n", lindex(n + 2, 0), c0, lindex(n, 0));
			fprintf(fp, "               vmulps         C%03i, r2, %%ymm9\n", constant(c0));
			fprintf(fp, "               vmulps         %i(Y), %%ymm9, %%ymm9\n", 32 * lindex(n, 0));
			fprintf(fp, "               vmovaps        %%ymm9, %i(Y)\n", 32 * lindex(n + 2, 0));
		}
		//tprint("Y[%i] = fma(TCAST(%.20e) * z, Y[%i], Y[%i]);\n", lindex(n + 1, 0), c1, lindex(n, 0), lindex(n + 1, 0));
		for (int m = 1; m <= n; m++) {
			const double c0 = -double(1) / (double((n + 2) * (n + 2)) - double(m * m));
			const double c1 = double(2 * n + 1) / (double((n + 1) * (n + 1)) - double(m * m));
			if (n + 2 <= P) {
				//tprint("ax0 = TCAST(%.20e) * r2;\n", c0);
				fprintf(fp, "               vmulps         C%03i, r2, ax0\n", constant(c0));
			}
			if (n + 2 <= P) {
				//	tprint("Y[%i] = ax0 * Y[%i];\n", lindex(n + 2, m), lindex(n, m));
				fprintf(fp, "               vmulps         %i(Y), ax0, %%ymm8\n", 32 * lindex(n, m));
				//	tprint("Y[%i] = ax0 * Y[%i];\n", lindex(n + 2, -m), lindex(n, -m));
				fprintf(fp, "               vmulps         %i(Y), ax0, %%ymm9\n", lindex(n, -m));
				fprintf(fp, "               vmovaps        %%ymm8, %i(Y)\n", 32 * lindex(n + 2, m));
				fprintf(fp, "               vmovaps        %%ymm9, %i(Y)\n", 32 * lindex(n + 2, -m));
			}
			if (n == m) {
				//	tprint("Y[%i] = z * Y[%i];\n", lindex(n + 1, m), lindex(n, m));
				fprintf(fp, "               vmulps         %i(Y), z, %%ymm8\n", 32 * lindex(n, m));
				//	tprint("Y[%i] = z * Y[%i];\n", lindex(n + 1, -m), lindex(n, -m));
				fprintf(fp, "               vmulps         %i(Y), z, %%ymm9\n", 32 * lindex(n, -m));
				fprintf(fp, "               vmovaps        %%ymm8, %i(Y)\n", 32 * lindex(n + 2, m));
				fprintf(fp, "               vmovaps        %%ymm9, %i(Y)\n", 32 * lindex(n + 2, -m));
			} else {
				//	tprint("ay0 = TCAST(%.20e) * z;\n", c1);
				fprintf(fp, "               vmulps         C%03i, z, ay0\n", constant(c1));
				//	tprint("Y[%i] = fma(ay0, Y[%i], Y[%i]);\n", lindex(n + 1, m), lindex(n, m), lindex(n + 1, m));
				fprintf(fp, "               vmovaps        %i(Y), %%ymm8\n", 32 * lindex(n + 1, m));
				fprintf(fp, "               vmovaps        %i(Y), %%ymm9\n", 32 * lindex(n + 1, -m));
				fprintf(fp, "               vfmadd231ps    %i(Y), ay0, %%ymm8\n", 32 * lindex(n, m));
				fprintf(fp, "               vfmadd231ps    %i(Y), ay0, %%ymm9\n", 32 * lindex(n, -m));
				fprintf(fp, "               vmovaps        %%ymm8, %i(Y)\n", 32 * lindex(n + 1, m));
				fprintf(fp, "               vmovaps        %%ymm9, %i(Y)\n", 32 * lindex(n + 1, -m));
				//	tprint("Y[%i] = fma(ay0, Y[%i], Y[%i]);\n", lindex(n + 1, -m), lindex(n, -m), lindex(n + 1, -m));
			}
		}
	}
	for (int i = 0; i < exp_sz(P); i += 16) {
		for (int j = 0; j < 16; j++) {
			int k = i + j;
			if (k < exp_sz(P)) {
				fprintf(fp, "               vmovaps        %i(M), %%ymm%i\n", 32 * k, j);
			}
		}
		for (int j = 0; j < 16; j++) {
			int k = i + j;
			if (k < exp_sz(P)) {
				fprintf(fp, "               vaddps         %i(Y), %%ymm%i, %%ymm%i\n", 32 * k, j, j);
			}
		}
		for (int j = 0; j < 16; j++) {
			int k = i + j;
			if (k < exp_sz(P)) {
				fprintf(fp, "               vmovaps        %%ymm%i, %i(M)\n", j, 32 * k);
			}
		}
	}
	fprintf(fp, "               mov            %%rsp, %%rbp\n");
	fprintf(fp, "               pop            %%rbp\n");
	fprintf(fp, "               ret\n");
	flush_constants(fp);
}

int greens_ewald(FILE* fp, int P) {
	int flops = 0;

	fprintf(fp, "#define        x0             %%ymm0\n");
	fprintf(fp, "#define        y0             %%ymm1\n");
	fprintf(fp, "#define        z0             %%ymm2\n");
	fprintf(fp, "#define        x              %%ymm3\n");
	fprintf(fp, "#define        y              %%ymm4\n");
	fprintf(fp, "#define        z              %%ymm5\n");
	fprintf(fp, "#define        r2             %%ymm6\n");
	fprintf(fp, "#define        r              %%ymm7\n");
	fprintf(fp, "#define        gamma1         %%ymm8\n");
	fprintf(fp, "#define        gamma          %%ymm9\n");
	fprintf(fp, "#define        sw             %%ymm10\n");
	fprintf(fp, "#define        xpow           %%ymm11\n");
	fprintf(fp, "#define        tmp            %%ymm12\n");
	fprintf(fp, "#define        tmp2           %%ymm13\n");
	fprintf(fp, "#define        hdotx          %%ymm14\n");
	fprintf(fp, "#define        cosphi         %%ymm15\n");
	fprintf(fp, "#define        sinphi         %%ymm8\n");
	fprintf(fp, "#define        G              %%rdi\n");
	fprintf(fp, "#define        Gr             %%rsi\n");
	fprintf(fp, "#define        erfc0          -8(%%rbp)\n");
	fprintf(fp, "#define        exp0           -16(%%rbp)\n");
	fprintf(fp, "#define        xfac           -24(%%rbp)\n");
	fprintf(fp, "\n");
	fprintf(fp, "               .global        greens_ewald\n");
	fprintf(fp, "\n");
	fprintf(fp, "               .text\n");
	fprintf(fp, "\n");
	fprintf(fp, "greens_ewald:  push           %%rbp\n");
	fprintf(fp, "               mov            %%rsp, %%rbp\n");
	fprintf(fp, "               sub            $%i, %%rsp\n", 32 * exp_sz(P) + 24);
	fprintf(fp, "               and            $0xffffffffffffffc0, %%rsp\n");
	fprintf(fp, "               mov            %%rsp, Gr\n");

	constexpr double alpha = 2.0;
	//tprint("T Gr[%i];\n", exp_sz(P));
	bool first = true;
	//tprint("T sw;\n");
	for (int ix = -3; ix <= 3; ix++) {
		for (int iy = -3; iy <= 3; iy++) {
			for (int iz = -3; iz <= 3; iz++) {
				if (ix * ix + iy * iy + iz * iz > 3.46 * 3.46) {
					continue;
				}
				//tprint("{\n");
				if (ix == 0) {
					//tprint("const T& x = x0;\n");
					fprintf(fp, "               vmovaps        x0, x\n");
				} else {
					//tprint("const T x = x0 - T(%i);\n", ix);
					fprintf(fp, "               vmovaps        x0, x\n");
					fprintf(fp, "               vsubps         C%03i, x, x\n", constant(double(ix)));
				}
				if (iy == 0) {
					//tprint("const T& y = y0;\n");
					fprintf(fp, "               vmovaps        y0, y\n");
				} else {
					//tprint("const T y = y0 - T(%i);\n", iy);
					fprintf(fp, "               vmovaps        y0, y\n");
					fprintf(fp, "               vsubps         C%03i, y, y\n", constant(double(iy)));
				}
				if (iz == 0) {
					//tprint("const T& z = z0;\n");
					fprintf(fp, "               vmovaps        z0, z\n");
				} else {
					//tprint("const T z = z0 - T(%i);\n", iz);
					fprintf(fp, "               vmovaps        z0, z\n");
					fprintf(fp, "               vsubps         C%03i, z, z\n", constant(double(iz)));
				}
				//tprint("const T r2 = x * x + y * y + z * z;\n");
				fprintf(fp, "               vmulps         z, z, r2\n");
				fprintf(fp, "               vfmadd231ps    y, y, r2\n");
				fprintf(fp, "               vfmadd231ps    x, x, r2\n");
				//tprint("const T r = sqrt(r2);\n");
				fprintf(fp, "               vsqrtps        r2, r\n");
				//tprint("greens_%s_P%i(Gr, x, y, z);\n", type, P);
				fprintf(fp, "               sub            $256, %%rsp\n");
				fprintf(fp, "               vmovaps        x0, (%%rsp)\n");
				fprintf(fp, "               vmovaps        y0, 32(%%rsp)\n");
				fprintf(fp, "               vmovaps        z0, 64(%%rsp)\n");
				fprintf(fp, "               vmovaps        x, 96(%%rsp)\n");
				fprintf(fp, "               vmovaps        y, 128(%%rsp)\n");
				fprintf(fp, "               vmovaps        z, 160(%%rsp)\n");
				fprintf(fp, "               vmovaps        r, 192(%%rsp)\n");
				fprintf(fp, "               vmovaps        r2, 224(%%rsp)\n");
				fprintf(fp, "               push           %%rdi\n");
				fprintf(fp, "               vmovaps        x, %%ymm0\n");
				fprintf(fp, "               vmovaps        y, %%ymm1\n");
				fprintf(fp, "               vmovaps        z, %%ymm2\n");
				fprintf(fp, "               mov            Gr, %%rdi\n");
				fprintf(fp, "               call           greens\n");
				fprintf(fp, "               vmovaps        r, %%ymm0\n");
				fprintf(fp, "               vmulps         C%03i, %%ymm0, %%ymm0\n", constant(alpha));
				fprintf(fp, "               vmovaps        %%ymm0, xpow\n");
				fprintf(fp, "               call           erfc_simd\n");
				fprintf(fp, "               vmovaps        %%ymm0, erfc0\n");
				fprintf(fp, "               vmovaps        r, %%ymm0\n");
				fprintf(fp, "               vmulps         C%03i, %%ymm0, %%ymm0\n", constant(alpha));
				fprintf(fp, "               vmulps         %%ymm0, %%ymm0, %%ymm0\n");
				fprintf(fp, "               vmovaps        %%ymm0, xfac\n");
				fprintf(fp, "               vmulps         C%03i, %%ymm0, %%ymm0\n", constant(-1));
				fprintf(fp, "               call           exp_simd\n");
				fprintf(fp, "               vmovaps        %%ymm0, exp0\n");
				fprintf(fp, "               pop            %%rdi\n");
				fprintf(fp, "               vmovaps        (%%rsp), x0\n");
				fprintf(fp, "               vmovaps        32(%%rsp), y0\n");
				fprintf(fp, "               vmovaps        64(%%rsp), z0\n");
				fprintf(fp, "               vmovaps        96(%%rsp), x\n");
				fprintf(fp, "               vmovaps        128(%%rsp), y\n");
				fprintf(fp, "               vmovaps        160(%%rsp), z\n");
				fprintf(fp, "               vmovaps        192(%%rsp), r\n");
				fprintf(fp, "               vmovaps        224(%%rsp), r2\n");
				fprintf(fp, "               add            $256, %%rsp\n");
				//tprint("T gamma1 = T(%.16e) * erfc(T(%.16e) * r);\n", sqrt(M_PI), alpha);
				fprintf(fp, "               vmovaps        erfc0, gamma1\n");
				fprintf(fp, "               vmulps         gamma1, r, gamma1\n");
				fprintf(fp, "               vmulps         C%03i, gamma1, gamma1\n", constant(sqrt(M_PI)));
				double gamma0inv = 1.0f / sqrt(M_PI);
				//tprint("T gamma;\n");
				if (ix * ix + iy * iy + iz * iz == 0) {
					//tprint("sw = (x0 * x0 + y0 * y0 + z0 * z0) > T(0);\n");
					fprintf(fp, "               vmulps         z0, z0, sw\n");
					fprintf(fp, "               vfmadd231ps    y0, y0, sw\n");
					fprintf(fp, "               vfmadd231ps    x0, x0, sw\n");
					fprintf(fp, "               vfmadd231ps    x0, x0, sw\n");
					fprintf(fp, "               vcmpps         $14, C%03i, sw, sw\n", constant(0.0));
					fprintf(fp, "               vpand          ABS, sw, sw\n");
					fprintf(fp, "               vcvtdq2ps      sw, sw\n");
				}
				for (int l = 0; l <= P; l++) {
					//tprint("gamma = gamma1 * T(%.16e);\n", gamma0inv);
					fprintf(fp, "               vmulps         C%03i, gamma1, gamma\n", constant(gamma0inv));
					if (ix * ix + iy * iy + iz * iz == 0) {
						for (int m = -l; m <= l; m++) {
							//tprint("G[%i] -= sw*(gamma - T(%.1e)) * Gr[%i];\n", index(l, m), nonepow<double>(l), index(l, m));
							fprintf(fp, "               vsubps         C%03i, gamma, tmp\n", constant(l % 2 == 0 ? 1 : -1));
							fprintf(fp, "               vmulps         %i(Gr), tmp, tmp\n", 32 * lindex(l, m));
							fprintf(fp, "               vmovaps        %i(G), tmp2\n", 32 * lindex(l, m));
							fprintf(fp, "               vfnmadd231ps   sw, tmp, tmp2\n");
							fprintf(fp, "               vmovaps        tmp2, %i(G)\n", 32 * lindex(l, m));
						}
						if (l == 0) {
							//tprint("G[%i] += (T(1) - sw)*T(%.16e);\n", index(0, 0), (2) * alpha / sqrt(M_PI));
							fprintf(fp, "               vmovaps        C%03i, tmp\n", constant(1));
							fprintf(fp, "               vmovaps        %i(G), tmp2\n", 32 * lindex(0, 0));
							fprintf(fp, "               vsubps         sw, tmp, tmp\n");
							fprintf(fp, "               vmulps         C%03i, tmp, tmp\n", constant((2) * alpha / sqrt(M_PI)));
							fprintf(fp, "               vaddps         tmp, tmp2, tmp2\n");
							fprintf(fp, "               vmovaps        tmp2, %i(G)\n", 32 * lindex(0, 0));
						}

					} else {
						if (first) {
							fprintf(fp, "               vmulps         C%03i, gamma, gamma\n", constant(-1));
							for (int m = -l; m <= l; m++) {
								//tprint("G[%i] = -gamma * Gr[%i];\n", index(l, m), index(l, m));
								fprintf(fp, "               vmulps         %i(Gr), gamma, tmp\n", 32 * lindex(l, m));
								fprintf(fp, "               vmovaps        tmp, %i(G)\n", 32 * lindex(l, m));
							}
							fprintf(fp, "               vmulps         C%03i, gamma, gamma\n", constant(-1));
						} else {
							for (int m = -l; m <= l; m++) {
								//tprint("G[%i] -= gamma * Gr[%i];\n", index(l, m), index(l, m));
								fprintf(fp, "               vmovaps        %i(G), tmp\n", 32 * lindex(l, m));
								fprintf(fp, "               vfnmadd231ps   %i(Gr), gamma, tmp\n", 32 * lindex(l, m));
								fprintf(fp, "               vmovaps        tmp, %i(G)\n", 32 * lindex(l, m));
							}
						}
					}
					gamma0inv *= 1.0 / -(l + 0.5);
					if (l != P) {
						//tprint("gamma1 = T(%.16e) * gamma1 + xpow * exp0;\n", l + 0.5);
						fprintf(fp, "               vmulps         exp0, xpow, tmp\n");
						fprintf(fp, "               vfmadd213ps    C%03i, tmp, gamma1\n", constant(l + 0.5));
						if (l != P - 1) {
							//tprint("xpow *= xfac;\n");
							fprintf(fp, "               vmulps         xfac, xpow, xpow\n");
						}
					}
				}
				//tprint("}\n");
				first = false;
			}
		}
	}
	//tprint("T cosphi, sinphi, hdotx, phi;\n");

	for (int hx = -2; hx <= 2; hx++) {
		for (int hy = -2; hy <= 2; hy++) {
			for (int hz = -2; hz <= 2; hz++) {
				const int h2 = hx * hx + hy * hy + hz * hz;
				if (h2 <= 8 && h2 > 0) {
					const double h = sqrt(h2);
					bool init = false;
					if (hx) {
						if (hx == 1) {
							//tprint("hdotx %s= x0;\n", init ? "+" : "");
							if (!init) {
								fprintf(fp, "               vmovaps        x0, hdotx\n");
							} else {
								fprintf(fp, "               vaddps         x0, hdotx, hdotx\n");
							}
							init = true;
						} else if (hx == -1) {
							if (init) {
								//tprint("hdotx -= x0;\n");
							} else {
								//tprint("hdotx = -x0;\n");
							}
							if (!init) {
								fprintf(fp, "               vmulps         C%03i, x0, hdotx\n", constant(-1));
							} else {
								fprintf(fp, "               vsubps         x0, hdotx, hdotx\n");
							}
							init = true;
						} else {
							//tprint("hdotx %s= T(%i) * x0;\n", init ? "+" : "", hx);
							if (!init) {
								fprintf(fp, "               vmulps         C%03i, x0, hdotx\n", constant(hx));
							} else {
								fprintf(fp, "               vfmadd231ps    C%03i, x0, hdotx\n", constant(hx));
							}
							init = true;
						}
					}
					if (hy) {
						if (hy == 1) {
							//tprint("hdotx %s= y0;\n", init ? "+" : "");
							if (!init) {
								fprintf(fp, "               vmovaps        y0, hdotx\n");
							} else {
								fprintf(fp, "               vaddps         y0, hdotx, hdotx\n");
							}
							init = true;
						} else if (hy == -1) {
							if (init) {
								//tprint("hdotx -= y0;\n");
							} else {
								//tprint("hdotx = -y0;\n");
							}
							if (!init) {
								fprintf(fp, "               vmulps         C%03i, y0, hdotx\n", constant(-1));
							} else {
								fprintf(fp, "               vsubps         y0, hdotx, hdotx\n");
							}
							init = true;
						} else {
							//tprint("hdotx %s= T(%i) * y0;\n", init ? "+" : "", hy);
							if (!init) {
								fprintf(fp, "               vmulps         C%03i, y0, hdotx\n", constant(hy));
							} else {
								fprintf(fp, "               vfmadd231ps    C%03i, y0, hdotx\n", constant(hy));
							}
							init = true;
						}
					}
					if (hz) {
						if (hz == 1) {
							//tprint("hdotx %s= z0;\n", init ? "+" : "");
							if (!init) {
								fprintf(fp, "               vmovaps        z0, hdotx\n");
							} else {
								fprintf(fp, "               vaddps         z0, hdotx, hdotx\n");
							}
							init = true;
						} else if (hz == -1) {
							if (init) {
								//tprint("hdotx -= z0;\n");
							} else {
								//tprint("hdotx = -z0;\n");
							}
							if (!init) {
								fprintf(fp, "               vmulps         C%03i, z0, hdotx\n", constant(-1));
							} else {
								fprintf(fp, "               vsubps         z0, hdotx, hdotx\n");
							}
							init = true;
						} else {
							//tprint("hdotx %s= T(%i) * z0;\n", init ? "+" : "", hz);
							if (!init) {
								fprintf(fp, "               vmulps         C%03i, z0, hdotx\n", constant(hz));
							} else {
								fprintf(fp, "               vfmadd231ps    C%03i, z0, hdotx\n", constant(hz));
							}
							init = true;
						}
					}

					const auto G0 = spherical_singular_harmonic(P, (double) hx, (double) hy, (double) hz);
					//tprint("phi = T(%.16e) * hdotx;\n", 2.0 * M_PI);
					fprintf(fp, "               sub            $128, %%rsp\n");
					fprintf(fp, "               vmovaps        x0, (%%rsp)\n");
					fprintf(fp, "               vmovaps        y0, 32(%%rsp)\n");
					fprintf(fp, "               vmovaps        z0, 64(%%rsp)\n");
					fprintf(fp, "               vmovaps        hdotx, 96(%%rsp)\n");
					fprintf(fp, "               vmulps         C%03i, hdotx, %%ymm0\n", constant(2.0 * M_PI));
					fprintf(fp, "               call           cos_simd\n");
					fprintf(fp, "               vmovaps        %%ymm0, cosphi\n");
					fprintf(fp, "               vmovaps        96(%%rsp), hdotx\n");
					fprintf(fp, "               vmulps         C%03i, hdotx, %%ymm0\n", constant(2.0 * M_PI));
					fprintf(fp, "               call           sin_simd\n");
					fprintf(fp, "               vmovaps        %%ymm0, sinphi\n");
					//tprint("SINCOS(phi, &sinphi, &cosphi);\n");
					fprintf(fp, "               vmovaps        (%%rsp), x0\n");
					fprintf(fp, "               vmovaps        32(%%rsp), y0\n");
					fprintf(fp, "               vmovaps        64(%%rsp), z0\n");
					fprintf(fp, "               add            $128, %%rsp\n");
					double gamma0inv = 1.0f / sqrt(M_PI);
					double hpow = 1.f / h;
					double pipow = 1.f / sqrt(M_PI);
					for (int l = 0; l <= P; l++) {
						for (int m = 0; m <= l; m++) {
							double c0 = gamma0inv * hpow * pipow * exp(-h * h * double(M_PI * M_PI) / (alpha * alpha));
							std::string ax = "cosphi";
							std::string ay = "sinphi";
							int xsgn = 1;
							int ysgn = 1;
							if (l % 4 == 3) {
								std::swap(ax, ay);
								ysgn = -1;
							} else if (l % 4 == 2) {
								ysgn = xsgn = -1;
							} else if (l % 4 == 1) {
								std::swap(ax, ay);
								xsgn = -1;
							}
							if (G0[cindex(l, m)].real() != (0)) {
								//tprint("G[%i] += T(%.16e) * %s;\n", index(l, m), -xsgn * c0 * G0[cindex(l, m)].real(), ax.c_str());
								fprintf(fp, "               vmovaps        %i(G), tmp\n", 32 * lindex(l, m));
								fprintf(fp, "               vfmadd231ps    C%03i, %s, tmp\n", constant(-xsgn * c0 * G0[cindex(l, m)].real()), ax.c_str());
								fprintf(fp, "               vmovaps        tmp, %i(G)\n", 32 * lindex(l, m));
								if (m != 0) {
									//tprint("G[%i] += T(%.16e) * %s;\n", index(l, -m), -ysgn * c0 * G0[cindex(l, m)].real(), ay.c_str());
									fprintf(fp, "               vmovaps        %i(G), tmp\n", 32 * lindex(l, -m));
									fprintf(fp, "               vfmadd231ps    C%03i, %s, tmp\n", constant(-ysgn * c0 * G0[cindex(l, m)].real()), ay.c_str());
									fprintf(fp, "               vmovaps        tmp, %i(G)\n", 32 * lindex(l, -m));
								}
							}
							if (G0[cindex(l, m)].imag() != (0)) {
								flops += 2;
								//tprint("G[%i] += T(%.16e) * %s;\n", index(l, m), ysgn * c0 * G0[cindex(l, m)].imag(), ay.c_str());
								fprintf(fp, "               vmovaps        %i(G), tmp\n", 32 * lindex(l, m));
								fprintf(fp, "               vfmadd231ps    C%03i, %s, tmp\n", constant(ysgn * c0 * G0[cindex(l, m)].real()), ay.c_str());
								fprintf(fp, "               vmovaps        tmp, %i(G)\n", 32 * lindex(l, m));
								if (m != 0) {
									flops += 2;
									//tprint("G[%i] += T(%.16e) * %s;\n", index(l, -m), -xsgn * c0 * G0[cindex(l, m)].imag(), ax.c_str());
									fprintf(fp, "               vmovaps        %i(G), tmp\n", 32 * lindex(l, -m));
									fprintf(fp, "               vfmadd231ps    C%03i, %s, tmp\n", constant(-xsgn * c0 * G0[cindex(l, m)].real()), ax.c_str());
									fprintf(fp, "               vmovaps        tmp, %i(G)\n", 32 * lindex(l, -m));
								}
							}
						}
						gamma0inv /= l + 0.5f;
						hpow *= h * h;
						pipow *= M_PI;

					}
				}
			}
		}
	}

	//tprint("G[%i] = T(%.16e);\n", (P + 1) * (P + 1), (4.0 * M_PI / 3.0));
	fprintf(fp, "               vmovaps        C%03i, %%ymm0\n", constant(4.0 * M_PI / 3.0));
	fprintf(fp, "               vmovaps        %%ymm0, %i(G)\n", (P + 1) * (P + 1));
	//tprint("/* flops = %i */\n", flops);
	//tprint("}");
	//tprint("\n");
	fprintf(fp, "               mov            %%rbp, %%rsp\n");
	fprintf(fp, "               pop            %%rbp\n");
	fprintf(fp, "               ret\n");
	flush_constants(fp);
	fprintf(fp, "ABS:           .long          0x7fffffff\n");
	fprintf(fp, "               .long          0x7fffffff\n");
	fprintf(fp, "               .long          0x7fffffff\n");
	fprintf(fp, "               .long          0x7fffffff\n");
	fprintf(fp, "               .long          0x7fffffff\n");
	fprintf(fp, "               .long          0x7fffffff\n");
	fprintf(fp, "               .long          0x7fffffff\n");
	fprintf(fp, "               .long          0x7fffffff\n");

	return flops;

}

struct entry_t {
	int l;
	int m;
	int o;
	int s;
	bool f;
};

const auto cmp = [](entry_t a, entry_t b) {
	if( a.l < b.l ) {
		return true;
	} else if( a.l > b.l ) {
		return false;
	} else {
		return a.m < b.m;
	}
};

void M2L(FILE* fp, int P) {
	fprintf(fp, "\n");
	fprintf(fp, "               .global        M2L\n");
	fprintf(fp, "\n");
	fprintf(fp, "               .text\n");
	fprintf(fp, "\n");
	fprintf(fp, "M2L:     \n");
	std::vector < std::queue < entry_t >> entries(exp_sz(P));
	bool first[(P + 1) * (P + 1)];
	for (int n = 0; n < (P + 1) * (P + 1); n++) {
		first[n] = true;
	}
	for (int n = 0; n <= P; n++) {
		for (int m = 0; m <= n; m++) {
			bool iload = true;
			bool rload = true;
			//tprint_new_chain();
			const int kmax = std::min(P - n, P - 1);
			for (int k = 0; k <= kmax; k++) {
				const int lmin = std::max(-k, -n - k - m);
				const int lmax = std::min(k, n + k - m);
				for (int l = lmin; l <= lmax; l++) {
					bool greal = false;
					bool mreal = false;
					int gxsgn = 1;
					int gysgn = 1;
					int mxsgn = 1;
					int mysgn = 1;
					int gxstr = -1;
					int gystr = -1;
					int mxstr = -1;
					int mystr = -1;
					if (m + l > 0) {
						gxstr = lindex(n + k, m + l);
						gystr = lindex(n + k, -m - l);
					} else if (m + l < 0) {
						if (abs(m + l) % 2 == 0) {
							gxstr = lindex(n + k, -m - l);
							gystr = lindex(n + k, m + l);
							gysgn = -1;
						} else {
							gxstr = lindex(n + k, -m - l);
							gystr = lindex(n + k, m + l);
							gxsgn = -1;
						}
					} else {
						greal = true;
						gxstr = lindex(n + k, 0);
					}
					if (l > 0) {
						mxstr = lindex(k, l);
						mystr = lindex(k, -l);
						mysgn = -1;
					} else if (l < 0) {
						if (l % 2 == 0) {
							mxstr = lindex(k, -l);
							mystr = lindex(k, l);
						} else {
							mxstr = lindex(k, -l);
							mystr = lindex(k, l);
							mxsgn = -1;
							mysgn = -1;
						}
					} else {
						mreal = true;
						mxstr = lindex(k, 0);
					}
					const auto csgn = [](int i) {
						return i > 0 ? '+' : '-';
					};
					const auto add_work = [n,&first, &gystr,&entries](int sgn, int m, int mstr, int gstr) {
						entry_t entry;
						entry.l = lindex(n, m);
						entry.m = mstr;
						entry.o = gstr;
						entry.s = sgn;
						entry.f = first[entry.l];
						first[entry.l] = false;
						entries[entry.l].push(entry);
					};
					if (!mreal) {
						if (!greal) {
							add_work(mxsgn * gxsgn, m, mxstr, gxstr);
							add_work(-mysgn * gysgn, m, mystr, gystr);
							if (m > 0) {
								add_work(mysgn * gxsgn, -m, mystr, gxstr);
								add_work(mxsgn * gysgn, -m, mxstr, gystr);
							}
						} else {
							add_work(mxsgn * gxsgn, m, mxstr, gxstr);
							if (m > 0) {
								add_work(mysgn * gxsgn, -m, mystr, gxstr);
							}
						}
					} else {
						if (!greal) {
							add_work(mxsgn * gxsgn, m, mxstr, gxstr);
							if (m > 0) {
								add_work(mxsgn * gysgn, -m, mxstr, gystr);
							}
						} else {
							add_work(mxsgn * gxsgn, m, mxstr, gxstr);
						}
					}
				}
			}
		}
	}
	int npar = 8;
	std::vector < std::queue < entry_t >> qs(npar);
	std::vector < std::queue < std::string >> cmds(npar);
	bool done = false;
	int iter = 0;
	while (!done) {
		done = true;
		char* str;
		for (int i = 0; i < npar; i++) {
			if (qs[i].size() == 0) {
				if (entries.size()) {
					done = false;
					qs[i] = std::move(entries.front());
					for (int l = 0; l < entries.size() - 1; l++) {
						entries[l] = std::move(entries[l + 1]);
					}
					entries.pop_back();
					asprintf(&str, "               vmovaps        %i(%%rdi), %%ymm%i\n", 32 * qs[i].front().l, 2 * i);
					cmds[i].push(str);
					free(str);
				}
			} else {
				done = false;
			}
			if (qs[i].size()) {
				auto& q = qs[i];
				auto e = q.front();
				asprintf(&str, "               vmovaps        %i(%%rdx), %%ymm%i\n", 32 * e.o, 2 * i + 1);
				cmds[i].push(str);
				free(str);
				if (e.f) {
					asprintf(&str, "               vmulps         %i(%%rsi), %%ymm%i, %%ymm%i\n", 32 * e.m, 2 * i + 1, 2 * i);
					cmds[i].push(str);
					free(str);
					if (e.s != 1) {
						abort();
					}
				} else {
					if (e.s == 1) {
						asprintf(&str, "               vfmadd231ps    %i(%%rsi), %%ymm%i, %%ymm%i\n", 32 * e.m, 2 * i + 1, 2 * i);
						cmds[i].push(str);
						free(str);
					} else {
						asprintf(&str, "               vfnmadd231ps   %i(%%rsi), %%ymm%i, %%ymm%i\n", 32 * e.m, 2 * i + 1, 2 * i);
						cmds[i].push(str);
						free(str);
					}

				}
				q.pop();
				if (!q.size()) {
					asprintf(&str, "               vmovaps        %%ymm%i, %i(%%rdi)\n", 2 * i, 32 * e.l);
					cmds[i].push(str);
					free(str);
				}
			}
		}
	}
	done = false;
	while (!done) {
		iter++;
		done = true;
		for (int i = 0; i < std::min(iter, npar); i++) {
			if (cmds[i].size()) {
				done = false;
				fprintf(fp, "%s", cmds[i].front().c_str());
				cmds[i].pop();
			}

		}
	}
	fprintf(fp, "               ret\n");
}

double factorial(int n) {
	return n == 0 ? 1.0 : n * factorial(n - 1);
}

void L2L_z(FILE* fp, int P, int Q) {
	const char* two = P == Q ? "" : "2";
	//tprint("Y[0] = %s;\n", var);
//	tprint("Y[1] = %s * %s;\n", var, var);
	fprintf(fp, "               vmovaps        %%ymm0, %%ymm15\n");
	fprintf(fp, "               vmulps         %%ymm0, %%ymm0, %%ymm14\n");
	int p = -1;
	for (int i = 2; i < P; i++) {
		const int n = (i - 1) / 2;
		const int m = (i - 1) - n;
//		tprint("Y[%i] = Y[%i] * Y[%i];\n", i, n, m);
		fprintf(fp, "               vmulps         %%ymm%i, %%ymm%i, %%ymm%i\n", 15 - m, 15 - n, 15 - i);
	}
	for (int i = 1; i < P; i++) {
//		tprint("Y[%i] *= TCAST(%0.20e);\n", i, 1.0 / factorial(i + 1));
		fprintf(fp, "               vmulps         C%03i, %%ymm%i, %%ymm%i\n", constant(1.0 / factorial(i + 1)), 15 - i, 15 - i);
	}
	for (int n = 0; n <= Q; n++) {
		std::vector < std::queue < entry_t >> cmds(exp_sz(Q));
		for (int m = -n; m <= n; m++) {
			for (int k = 1; k <= P - n; k++) {
				if (abs(m) > n + k) {
					continue;
				}
				entry_t e;
				e.m = lindex(n + k, m);
				e.l = lindex(n, m);
				e.o = k - 1;
				cmds[e.l].push(e);
				//tprint("L%s[%i] = fma(Y[%i], L[%i], L%s[%i]);\n", two, lindex(n, m), k - 1, lindex(n + k, m), two, lindex(n, m));
			}
		}
		int npar = (15 - P) / 2;
		std::vector < std::queue < entry_t >> qs(npar);
		bool done = false;
		while (!done) {
			done = true;
			for (int m = 0; m < npar; m++) {
				while (qs[m].size() == 0 && cmds.size()) {
					qs[m] = std::move(cmds.back());
					cmds.pop_back();
					if (qs[m].size()) {
						fprintf(fp, "               vmovaps        %i(%%rsi), %%ymm%i\n", 32 * qs[m].front().l, 2 * m);
						done = false;
					}
				}
			}
			for (int m = 0; m < npar; m++) {
				if (!qs[m].empty()) {
					done = false;
					auto e = qs[m].front();
					fprintf(fp, "               vmovaps        %i(%%rsi), %%ymm%i\n", 32 * e.m, 2 * m + 1);
				}
			}
			for (int m = 0; m < npar; m++) {
				if (!qs[m].empty()) {
					auto& e = qs[m].front();
					fprintf(fp, "               vfmadd231ps    %%ymm%i, %%ymm%i, %%ymm%i\n", 15 - e.o, 2 * m + 1, 2 * m);
				}
			}
			for (int m = 0; m < npar; m++) {
				int l;
				if (!qs[m].empty()) {
					l = qs[m].front().l;
					qs[m].pop();
					if (qs[m].empty()) {
						fprintf(fp, "               vmovaps        %%ymm%i, %i(%%rdi)\n", 2 * m, 32 * l);
					}
				}
			}
		}
	}
}

void z_rot2(FILE* fp, int P, bool stage) {
	std::pair<int, int> range;
	range.first = 1;
	range.second = P;

	if (stage) {
		//tprint("r0[0] = cosphi;\n");
		fprintf(fp, "               vmovaps        %%ymm0, (%%r8)\n");
		//tprint("ip[0] = sinphi;\n");
		fprintf(fp, "               vmovaps        %%ymm1, (%%r9)\n");
		for (int m = 1; m < P; m++) {
			if (m == 1) {
				fprintf(fp, "               vmovaps        %%ymm0, %%ymm4\n");
				fprintf(fp, "               vmovaps        %%ymm1, %%ymm5\n");
			} else {
				fprintf(fp, "               vmovaps        %%ymm2, %%ymm4\n");
				fprintf(fp, "               vmovaps        %%ymm3, %%ymm5\n");
			}
			//tprint("r0[%i] = in[%i] * sinphi;\n", m, m - 1);
			fprintf(fp, "               vmulps         %%ymm5, %%ymm1, %%ymm2\n");
			//tprint("ip[%i] = ip[%i] * cosphi;\n", m, m - 1);
			fprintf(fp, "               vmulps         %%ymm5, %%ymm0, %%ymm3\n");
			//tprint("r0[%i] = fma(r0[%i], cosphi, r0[%i]);\n", m, m - 1, m);
			fprintf(fp, "               vfnmadd231ps   %%ymm4, %%ymm0, %%ymm2\n");
			//tprint("ip[%i] = fma(r0[%i], sinphi, ip[%i]);\n", m, m - 1, m);
			fprintf(fp, "               vfmadd231ps    %%ymm4, %%ymm1, %%ymm3\n");
			fprintf(fp, "               vmovaps        %%ymm2, %i(%%r8)\n", 32 * m);
			fprintf(fp, "               vmovaps        %%ymm3, %i(%%r9)\n", 32 * m);
		}
	}
	int mmin = 1;
	bool initR = true;
	std::vector<int> set_nan;
	using cmd_t = std::pair<int,std::string>;
	std::vector<cmd_t> cmds;
//	cmds.push_back(std::make_pair(0, print2str("%s[0] = %s[0];\n", dst, src)));
	for (int m = range.first; m <= range.second; m++) {
		//	cmds.push_back(std::make_pair(index(m, 0), print2str("%s[%i] = %s[%i];\n", dst, index(m, 0), src, index(m, 0))));
		fprintf(fp, "               vmovaps        %i(%%r8), %%ymm14\n", 32 * (m - 1));
		fprintf(fp, "               vmovaps        %i(%%r9), %%ymm15\n", 32 * (m - 1));
		for (int l = m; l <= P; l++) {
			bool write_ronly = false;
			if (stage) {
				write_ronly = m % 2 != l % 2;
			}
			bool sw = false;
			if (write_ronly) {
				//	cmds.push_back(std::make_pair(index(l, -m), print2str("%s[%i] = %s[%i] * in[%i];\n", dst, index(l, m), src, index(l, -m), m - 1)));
				//	cmds.push_back(std::make_pair(index(l, m), print2str("%s[%i] = fma(%s[%i], r0[%i], %s[%i]);\n", dst, index(l, m), src, index(l, m), m - 1, dst, index(l, m))));
				fprintf(fp, "               vmulps         %i(%%rsi), %%ymm15, %%ymm1\n", 32 * lindex(l, -m));
				if (stage) {
					fprintf(fp, "               vfmsub231ps    %i(%%rsi), %%ymm14, %%ymm1\n", 32 * lindex(l, m));
				} else {
					fprintf(fp, "               vfmadd231ps    %i(%%rsi), %%ymm14, %%ymm1\n", 32 * lindex(l, m));
				}
				fprintf(fp, "               vmovaps        %%ymm1, %i(%%rdi)\n", 32 * lindex(l, m));
			} else {
				fprintf(fp, "               vmovaps        %i(%%rsi), %%ymm2\n", 32 * lindex(l, m));
				fprintf(fp, "               vmovaps        %i(%%rsi), %%ymm3\n", 32 * lindex(l, -m));

				fprintf(fp, "               vmulps         %%ymm3, %%ymm14, %%ymm0\n");
				fprintf(fp, "               vmulps         %%ymm2, %%ymm15, %%ymm1\n");
				if (stage) {
					fprintf(fp, "               vfmsub231ps    %%ymm2, %%ymm14, %%ymm0\n");
					fprintf(fp, "               vfmadd231ps    %%ymm3, %%ymm15, %%ymm1\n");
				} else {
					fprintf(fp, "               vfmadd231ps    %%ymm2, %%ymm14, %%ymm0\n");
					fprintf(fp, "               vfnmadd231ps   %%ymm3, %%ymm15, %%ymm1\n");
				}
				fprintf(fp, "               vmovaps        %%ymm0, %i(%%rdi)\n", 32 * lindex(l, m));
				fprintf(fp, "               vmovaps        %%ymm1, %i(%%rdi)\n", 32 * lindex(l, -m));
				//	cmds.push_back(std::make_pair(index(l, -m), print2str("%s[%i] = %s[%i] * in[%i];\n", dst, index(l, m), src, index(l, -m), m - 1)));
				//	cmds.push_back(std::make_pair(index(l, -m), print2str("%s[%i] = %s[%i] * r0[%i];\n", dst, index(l, -m), src, index(l, -m), m - 1)));
				//	cmds.push_back(std::make_pair(index(l, m), print2str("%s[%i] = fma(%s[%i], r0[%i], %s[%i]);\n", dst, index(l, m), src, index(l, m), m - 1, dst, index(l, m))));
				//	cmds.push_back(std::make_pair(index(l, m), print2str("%s[%i] = fma(%s[%i], ip[%i], %s[%i]);\n", dst, index(l, -m), src, index(l, m), m - 1, dst, index(l, -m))));
			}

		}
	}
	std::sort(cmds.begin(), cmds.end(), [](const cmd_t& a, const cmd_t& b) {
		return a.first < b.first;
	});
	for (auto cmd : cmds) {
		//tprint("%s", cmd.second.c_str());
	}
}

//xz_swap2(P, "L", "L2", true, PRE2, "L2P");
template<class T>
T nonepow(int m) {
	return m % 2 == 0 ? double(1) : double(-1);
}

void xz_swap2(FILE* fp, int P, bool stage) {

	auto index = lindex;
	struct entry_t {
		int s;
		int d;
		bool first;
		double c0;
	};
	std::vector < std::queue < entry_t >> ops(exp_sz(P));
	for (int n = 1; n <= P; n++) {
		int lmax = n;
		int mmax = n;
		if (stage) {
			mmax = 1;
		}
		for (int m = 0; m <= n; m++) {
			if (m <= mmax) {
				for (int l = 0; l <= lmax; l++) {
					bool flag = false;
					double r = l == 0 ? brot(n, m, 0) : brot(n, m, l) + nonepow<double>(l) * brot(n, m, -l);
					double i = l == 0 ? 0.0 : brot(n, m, l) - nonepow<double>(l) * brot(n, m, -l);
					if (r != 0.0) {
						entry_t e;
						e.d = index(n, m);
						e.s = index(n, l);
						e.c0 = r;
						e.first = ops[n + m].empty();
						ops[index(n, m)].push(e);
					}
					if (i != 0.0 && m != 0) {
						entry_t e;
						e.d = index(n, n - 2 * m);
						e.c0 = i;
						e.s = index(n, -l);
						e.first = ops[n - m].empty();
						ops[index(n, n - 2 * m)].push(e);
					}
				}
			}
		}
	}
	for (int l = 0; l < exp_sz(P); l++) {
		if (ops[l].size()) {
			entry_t e;
			e.d = ops[l].back().d;
			ops[l].push(e);
		}
	}
	int npar = 8;
	while (ops.size()) {
		int sz = std::min((int) ops.size(), npar);
		for (int l = 0; l < sz; l++) {
			while (ops.size() > l && ops[l].size() == 0) {
				ops[l] = std::move(ops.back());
				ops.pop_back();
			}
		}
		for (int l = 0; l < sz; l++) {
			if (!ops[l].size()) {
				continue;
			}
			auto& op = ops[l].front();
			bool st = ops[l].size() == 1;
			if (st) {
				fprintf(fp, "               vmovaps        %%ymm%i, %i(%%rdi)\n", 2 * l, 32 * op.d);
			} else {
				if (close2(1.0, op.c0)) {
					if (op.first) {
						fprintf(fp, "               vmovaps        %i(%%rsi), %%ymm%i\n", 32 * op.s, 2 * l);
					} else {
						fprintf(fp, "               vaddps         %i(%%rsi), %%ymm%i, %%ymm%i\n", 32 * op.s, 2 * l, 2 * l);
					}
				} else {
					fprintf(fp, "               vmovaps        %i(%%rsi), %%ymm%i\n", 32 * op.s, 2 * l + 1);
					if (op.first) {
						fprintf(fp, "               vmulps         C%03i, %%ymm%i, %%ymm%i\n", constant(op.c0), 2 * l + 1, 2 * l);
					} else {
						fprintf(fp, "               vfmadd231ps    C%03i, %%ymm%i, %%ymm%i\n", constant(op.c0), 2 * l + 1, 2 * l);
					}
				}
			}
			ops[l].pop();
		}
	}

}

void L2P(FILE* fp, int P) {
	fprintf(fp, "#define        x              -32(%%rbp)\n");
	fprintf(fp, "#define        y              -64(%%rbp)\n");
	fprintf(fp, "#define        z              -96(%%rbp)\n");
	fprintf(fp, "#define        cosphi         -128(%%rbp)\n");
	fprintf(fp, "#define        sinphi         -160(%%rbp)\n");
	fprintf(fp, "#define        R              -192(%%rbp)\n");
	fprintf(fp, "#define        Lsrc           -224(%%rbp)\n");
	fprintf(fp, "#define        Ldst           -232(%%rbp)\n");
	fprintf(fp, "#define        Lwrk           -240(%%rbp)\n");
	fprintf(fp, "\n");
	fprintf(fp, "               .global        L2P\n");
	fprintf(fp, "\n");
	fprintf(fp, "               .text\n");
	fprintf(fp, "\n");
	fprintf(fp, "L2P:           push           %%rbp\n");
	fprintf(fp, "               mov            %%rsp, %%rbp\n");
	fprintf(fp, "               sub            $%i, %%rsp\n", 248 + exp_sz(P) * 32);
	fprintf(fp, "               and            $0xffffffffffffffc0, %%rsp\n");
	fprintf(fp, "               mov            %%rsp, Lwrk\n");
	fprintf(fp, "               sub            $%i, %%rsp\n", 32 * P);
	fprintf(fp, "               mov            %%rsp, %%r8\n");
	fprintf(fp, "               sub            $%i, %%rsp\n", 32 * P);
	fprintf(fp, "               mov            %%rsp, %%r9\n");
	fprintf(fp, "               mov            %%rsi, Lsrc\n");
	fprintf(fp, "               mov            %%rdi, Ldst\n");
	fprintf(fp, "               vmovups        %%ymm0, x\n");
	fprintf(fp, "               vmovups        %%ymm1, y\n");
	fprintf(fp, "               vmovups        %%ymm2, z\n");

//	tprint("Rzero = TCONVERT( R2 < TCAST(%.20e) );\n", tiny());
	fprintf(fp, "               vmulps         %%ymm0, %%ymm0, %%ymm3\n");
	fprintf(fp, "               vfmadd231ps    %%ymm1, %%ymm1, %%ymm3\n");
	fprintf(fp, "               vfmadd231ps    %%ymm2, %%ymm2, %%ymm3\n");
	fprintf(fp, "               vcmpps         $1, C%03i, %%ymm3, %%ymm4\n", constant(1e-6));
	fprintf(fp, "               vpand          ABS, %%ymm4, %%ymm4\n");
//	tprint("tmp1 = R2 + Rzero;\n");
	fprintf(fp, "               vaddps         %%ymm3, %%ymm4, %%ymm5\n");
//	tprint("Rinv = rsqrt(tmp1);\n");
	fprintf(fp, "               vrsqrtps       %%ymm5, %%ymm6\n");
//	tprint("R = (TCAST(1) - Rzero) / Rinv;\n");
	fprintf(fp, "               vmovaps        C%03i, %%ymm7\n", constant(1));
	fprintf(fp, "               vsubps         %%ymm4, %%ymm7, %%ymm8\n");
	fprintf(fp, "               vdivps         %%ymm6, %%ymm8, %%ymm9\n");
	fprintf(fp, "               vmovaps        %%ymm9, R\n");
//	tprint("cosphi = fma(x, Rinv, Rzero);\n");
	fprintf(fp, "               vfmadd231ps    x, %%ymm6, %%ymm4\n");
//	tprint("sinphi = -y * Rinv;\n");*/
	fprintf(fp, "               vmulps         y, %%ymm6, %%ymm5\n");
	fprintf(fp, "               vmulps         C%03i, %%ymm5, %%ymm5\n", constant(-1.0));
	fprintf(fp, "               vmovups        %%ymm4, cosphi\n");
	fprintf(fp, "               vmovups        %%ymm5, sinphi\n");
	fprintf(fp, "               vmovups        z, %%ymm0\n");
	fprintf(fp, "               mov            Lsrc, %%rsi\n");
	fprintf(fp, "               mov            Lsrc, %%rdi\n");
	L2L_z(fp, P, P);
	fprintf(fp, "               mov            Lwrk, %%rsi\n");
	fprintf(fp, "               mov            Ldst, %%rdi\n");
	z_rot2(fp, P, true);
	fprintf(fp, "               vmovups        cosphi, %%ymm0\n");
	fprintf(fp, "               vmovups        sinphi, %%ymm1\n");
	fprintf(fp, "               mov            Lwrk, %%rdi\n");
	fprintf(fp, "               mov            Ldst, %%rsi\n");
	xz_swap2(fp, P, true);
	fprintf(fp, "               mov            Ldst, %%rdi\n");
	fprintf(fp, "               mov            Lwrk, %%rsi\n");
	fprintf(fp, "               vmovaps        (%%rsi), %%ymm0\n");
	fprintf(fp, "               vmovaps        32(%%rsi), %%ymm1\n");
	fprintf(fp, "               vmovaps        64(%%rsi), %%ymm2\n");
	fprintf(fp, "               vmovaps        128(%%rsi), %%ymm3\n");
	fprintf(fp, "               vmovaps        %%ymm0, (%%rdi)\n");
	fprintf(fp, "               vmovaps        %%ymm1, 32(%%rdi)\n");
	fprintf(fp, "               vmovaps        %%ymm2, 64(%%rdi)\n");
	fprintf(fp, "               vmovaps        %%ymm3, 128(%%rdi)\n");
	fprintf(fp, "               vmovups        R, %%ymm0\n");
	L2L_z(fp, P, 1);
	fprintf(fp, "               vmovups        cosphi, %%ymm0\n");
	fprintf(fp, "               vmovups        sinphi, %%ymm1\n");
	fprintf(fp, "               mov            Lwrk, %%rdi\n");
	fprintf(fp, "               mov            Ldst, %%rsi\n");
	xz_swap2(fp, 1, false);
	fprintf(fp, "               mov            Lwrk, %%rsi\n");
	fprintf(fp, "               mov            Ldst, %%rdi\n");
	z_rot2(fp, 1, false);
	fprintf(fp, "               mov            %%rbp, %%rsp\n");
	fprintf(fp, "               pop            %%rbp\n");
	fprintf(fp, "               ret\n");
	flush_constants(fp);
	fprintf(fp, "ABS:           .long          0x7fffffff\n");
	fprintf(fp, "               .long          0x7fffffff\n");
	fprintf(fp, "               .long          0x7fffffff\n");
	fprintf(fp, "               .long          0x7fffffff\n");
	fprintf(fp, "               .long          0x7fffffff\n");
	fprintf(fp, "               .long          0x7fffffff\n");
	fprintf(fp, "               .long          0x7fffffff\n");
	fprintf(fp, "               .long          0x7fffffff\n");
}

int main() {
	FILE* fp = fopen("greens.S", "wt");
	greens(fp, ORDER);
	fclose(fp);
	fp = fopen("greens_ewald.S", "wt");
	greens_ewald(fp, ORDER);
	fclose(fp);
	fp = fopen("M2L.S", "wt");
	M2L(fp, ORDER);
	fclose(fp);
	fp = fopen("P2M.S", "wt");
	P2M(fp, ORDER - 1);
	fclose(fp);
	fp = fopen("L2P.S", "wt");
	L2P(fp, ORDER);
	fclose(fp);

}
