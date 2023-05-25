#undef NDEBUG
#include "assert.h"
#include "verifier.h"
#include "inputCircuit.hpp"

using namespace std;
using namespace virgo;

void test_field_arithmetic() {
	virgo::fieldElement a, b, c, d, e;
	for (int i = 0; i < 1000; i++) {
		a = fieldElement::random();
		b = fieldElement::random();
		d = 0;

		c = fieldElement::zero();
		assert(c == 0);

		c = fieldElement::one();
		assert(c == 1);

		c = a * b;
		d = b * a;
		assert(c == d);

		c = fieldElement::random();
		d = (a * b) * c;
		e = a * (b * c);
		assert(d == e);

		c = a.inv();
		d = a * c;
		assert(d == 1);

		c = -1;
		assert(c.isNegative());
		c = c + 1;
		assert(!c.isNegative());

		c = a * a;
		d = a.sqr();
		assert(c == d);

		for (int i = 1; i < fieldElement::maxOrder(); i++) {
			c = fieldElement::getRootOfUnity(i);
			for (int j = 1; j < i; j++) {
				c = c.sqr();
			}
			assert(c == -1);
		}
	}
	cout << "Field tests pass. " << endl;
	return;

	fieldElement::self_speed_test_add(100);
	fieldElement::self_speed_test_mult(100);
}

void DAG_to_layered(layeredCircuit &c, vector<DAG_gate *> in_circuit_dag) {
	const int repeat = 1;
    vector<u64> in_deg(in_circuit_dag.size());          // in degree
    vector<int> lyr_id(in_circuit_dag.size());          // the layer of each gate
    vector<u64> id_in_lyr(in_circuit_dag.size());       // the corresponding id within the layer
    vector<vector<u64>> edges(in_circuit_dag.size());   // the edges in the DAG

    // Topologically sorting
    queue<u64> q;
    for (u64 i = 0; i < in_circuit_dag.size(); ++i) {
        auto &g = *in_circuit_dag[i];
        if (g.input0.first == 'V') {
            ++in_deg[i];
            edges[g.input0.second].push_back(i);
        }
        if (g.input1.first == 'V') {
            ++in_deg[i];
            edges[g.input1.second].push_back(i);
        }
        if (g.ty == Input) {
            lyr_id[i] = 0;
            q.push(i);
        }
    }

    int max_lyr_id = 0;
    while (!q.empty()) {
        u64 u = q.front();
        q.pop();
        max_lyr_id = max(lyr_id[u], max_lyr_id);
        for (auto v: edges[u])
            if (!(--in_deg[v])) {
                q.push(v);
                lyr_id[v] = max(lyr_id[v], lyr_id[u] + 1);
            }
    }

    // build the circuit
    c.circuit.resize(max_lyr_id + 1);
    c.size = max_lyr_id + 1;

    for (u64 i = 0; i < in_circuit_dag.size(); ++i)
        id_in_lyr[i] = c.circuit[lyr_id[i]].size++;

    for (int i = 0; i < c.size; ++i)
        c.circuit[i].gates.resize(c.circuit[i].size);

    for (u64 i = 0; i < in_circuit_dag.size(); ++i) {
        int lg = lyr_id[i];
        u64 gid = id_in_lyr[i];
        auto &g = *in_circuit_dag[i];
        auto ty = g.ty, nty = ty;
        u64 in0 = g.input0.second;
        u64 in1 = g.input1.second;
        bool is_assert = g.is_assert;
        u64 u, v;
        F cst;

        switch (ty) {
            case Mul: case Add: case Xor:
                u = id_in_lyr[in0];
                v = id_in_lyr[in1];
                if (lyr_id[in0] < lg - 1) swap(u, v), swap(in0, in1);
                c.circuit[lg].gates[gid] = gate(ty, lyr_id[in1], u, v, F_ZERO, is_assert);
            break;
            case Sub:
                u = id_in_lyr[in0];
                v = id_in_lyr[in1];
                if (lyr_id[in0] < lg - 1) {
                    nty = AntiSub;
                    swap(u, v);
                    swap(in0, in1);
                }
                c.circuit[lg].gates[gid] = gate(nty, lyr_id[in1], u, v, F_ZERO, is_assert);
            break;
            case Naab:
                u = id_in_lyr[in0];
                v = id_in_lyr[in1];
                if (lyr_id[in0] < lg - 1) {
                    nty = AntiNaab;
                    swap(u, v);
                    swap(in0, in1);
                }
                c.circuit[lg].gates[gid] = gate(nty, lyr_id[in1], u, v, F_ZERO, is_assert);
            break;
            case Mulc: case Addc:
                u = id_in_lyr[in0];
                cst = F(in1);
                c.circuit[lg].gates[gid] = gate(ty, -1, u, 0, cst, is_assert);
            break;
            case Not: case Copy:
                u = id_in_lyr[in0];
                cst = F(in1);
                c.circuit[lg].gates[gid] = gate(ty, -1, u, 0, cst, is_assert);
            case Input:
                u = in0;
                c.circuit[lg].gates[gid] = gate(ty, -1, u, 0, F_ZERO, is_assert);
        }
    }

    // repeat the layer except the input for ${repeat} times
    for (int i = 1; i < c.size; ++i) {
        for (int j = 1; j < repeat; ++j)
            for (u64 k = 0; k < c.circuit[i].size; ++k) {
                auto &g = c.circuit[i].gates[k];
                c.circuit[i].gates.push_back(c.circuit[i].gates[k]);
                if (g.ty != Input && i > 1) g.u += j * c.circuit[i].size;
                if (g.ty == Add ||
                    g.ty == Mul ||
                    g.ty == Xor ||
                    g.ty == Sub ||
                    g.ty == AntiSub ||
                    g.ty == Naab ||
                    g.ty == AntiNaab)
                    g.v += j * c.circuit[i].size;
            }
        c.circuit[i].size *= repeat;
    }

    for (int i = 0; i <= max_lyr_id; ++i) {
        c.circuit[i].bitLength = (int) log2(c.circuit[i].size);
        if ((1ULL << c.circuit[i].bitLength) < c.circuit[i].size) ++c.circuit[i].bitLength;
    }
}

regex add_gate("P V[0-9]+ = V[0-9]+ \\+ V[0-9]+ E");
regex mult_gate("P V[0-9]+ = V[0-9]+ \\* V[0-9]+ E");
regex input_gate("P V[0-9]+ = I[0-9]+ E");
regex output_gate("P O[0-9]+ = V[0-9]+ E");
regex xor_gate("P V[0-9]+ = V[0-9]+ XOR V[0-9]+ E");
regex minus_gate("P V[0-9]+ = V[0-9]+ minus V[0-9]+ E");
regex naab_gate("P V[0-9]+ = V[0-9]+ NAAB V[0-9]+ E");
regex not_gate("P V[0-9]+ = V[0-9]+ NOT V[0-9]+ E");

smatch base_match;

DAG_gate *buildGate(vector<DAG_gate *> &in_circuit_dag, gateType ty, u64 tgt,
					u64 src0, u64 src1, bool has_constant) {
//	fprintf(stderr, "buildGate: tgt: %d, src0: %d, src1: %d, has_const: %d\n", tgt, src0, src1, (int) has_constant);
    DAG_gate *g = new DAG_gate();
    g->is_assert = false;
    g->ty = ty;
    g->input0 = make_pair((int)'V', src0);
    g->input1 = make_pair(has_constant ? (int)'S' : 'V', src1);
    if (tgt >= in_circuit_dag.size()) in_circuit_dag.resize(tgt + 1, nullptr);
    in_circuit_dag[tgt] = g;
    return g;
}

DAG_gate *buildInput(vector<DAG_gate *> &in_circuit_dag, u64 tgt, u64 src0) {
//	fprintf(stderr, "buildInput: tgt: %d, src0: %d\n", tgt, src0);
    DAG_gate *g = new DAG_gate();
    g->is_assert = false;
    g->ty = Input;
    g->input0 = make_pair((int)'S', src0);
    g->input1 = make_pair((int)'N', 0);
    if (tgt >= in_circuit_dag.size()) in_circuit_dag.resize(tgt + 1, nullptr);
    in_circuit_dag[tgt] = g;
    return g;
}

void parse(vector<DAG_gate *> &in_circuit_dag, ifstream &circuit_in) {
    string source_line;
    i64 tgt, src0, src1;
    while (getline(circuit_in, source_line)) {
        if (std::regex_match(source_line, base_match, add_gate)) {
            sscanf(source_line.c_str(), "P V%lld = V%lld + V%lld E", &tgt, &src0, &src1);
            buildGate(in_circuit_dag, Add, tgt, src0, src1, false);
        } else if (std::regex_match(source_line, base_match, mult_gate)) {
            sscanf(source_line.c_str(), "P V%lld = V%lld * V%lld E", &tgt, &src0, &src1);
            buildGate(in_circuit_dag, Mul, tgt, src0, src1, false);
        } else if (std::regex_match(source_line, base_match, input_gate)) {
            sscanf(source_line.c_str(), "P V%lld = I%lld E", &tgt, &src0);
            buildInput(in_circuit_dag, tgt, src0);
        } else if (std::regex_match(source_line, base_match, output_gate)) {
            sscanf(source_line.c_str(), "P O%lld = V%lld E", &tgt, &src0);
        } else if (std::regex_match(source_line, base_match, xor_gate)) {
            sscanf(source_line.c_str(), "P V%lld = V%lld XOR V%lld E", &tgt, &src0, &src1);
            buildGate(in_circuit_dag, Xor, tgt, src0, src1, false);
        } else if (std::regex_match(source_line, base_match, naab_gate)) {
            sscanf(source_line.c_str(), "P V%lld = V%lld NAAB V%lld E", &tgt, &src0, &src1);
            buildGate(in_circuit_dag, Naab, tgt, src0, src1, false);
        } else if (std::regex_match(source_line, base_match, minus_gate)) {
            sscanf(source_line.c_str(), "P V%lld = V%lld minus V%lld E", &tgt, &src0, &src1);
            buildGate(in_circuit_dag, Sub, tgt, src0, src1, false);
        } else if (std::regex_match(source_line, base_match, not_gate)) {
            sscanf(source_line.c_str(), "P V%lld = V%lld NOT V%lld E", &tgt, &src0, &src1);
            buildGate(in_circuit_dag, Not, tgt, src0, 0, true);
        } else {
            assert(false);
        }
    }
}

#define REL 1
#define PRIME 0x1ffffc0000001LL
#define ROOT  416204888522856

#define PRIME2 0x2a74200000001LL
#define ROOT2  186427948752465

F* public_array_prepare_generic(F *public_array, int log_length)
{
	F *q_coef_arr = new F[1 << log_length];
	int coef_slice_size = (1 << (log_length - log_slice_number));
	for(int i = 0; i < (1 << log_slice_number); ++i)
	{
		inverse_fast_fourier_transform(&public_array[i * coef_slice_size], coef_slice_size, coef_slice_size, F::getRootOfUnity(log_length - log_slice_number), &q_coef_arr[i * coef_slice_size]);
	}
	return q_coef_arr;
}

int main(int argc, char **argv) {
	layeredCircuit c;
	vector<DAG_gate *> in_circuit_dag;

	// Configure prime field.
	F::init(PRIME, ROOT);
	test_field_arithmetic();

    ifstream circuit_in(argv[REL]);
    parse(in_circuit_dag, circuit_in);
    DAG_to_layered(c, in_circuit_dag);

    fclose(stdin);

	c.subsetInit();
    prover p(c);
    verifier v(&p, c);
    v.verify();

    fprintf(stdout, "mult counter %d, add counter %d\n", F::multCounter, F::addCounter);

    // Configure another prime field.
    F::init(PRIME2, ROOT2);
    test_field_arithmetic();

    //prover p2(c);
    //verifier v2(&p2, c);
	
	//F inner_product_sum;
    //vector<F> processed, all_sum(slice_number + 1), output(1ULL << c.circuit[0].bitLength);
    //v2.public_array_prepare_generic(processed, output, c.circuit[0].bitLength);

	/* Prover. */
	//auto mask2 = vector<F>(1, F_ZERO);
	//auto merkle_root_l = p2.commit_private();
    //auto merkle_root_h = p2.commit_public(output, inner_product_sum, mask2, all_sum);

	/* Verifier. */
    //v2.verify(all_sum, processed, mask2, merkle_root_l, merkle_root_h);
    //fprintf(stdout, "mult counter %d, add counter %d\n", F::multCounter, F::addCounter);
	
	/* Now just the polynomial commitment. */
	poly_commit::poly_commit_prover prover;
	poly_commit::poly_commit_verifier verifier;

	verifier.p = &prover;

	F *all_sum = new F[slice_number + 1];
	int log_length = 8;
	auto all_pri_mask = vector<F>(1, F_ZERO);
	auto all_pub_mask = vector<F>(1, F_ZERO);

	/* Prover. */
	vector<F> public_array(1 << log_length), private_array(1 << log_length);
	for (int i = 0; i < (1 << log_length); i++) {
		private_array[i] = fieldElement::random();
		public_array[i] = fieldElement::random();
	}
	auto merkle_root_l = prover.commit_private_array(private_array.data(), log_length, all_pri_mask);
	auto inner_product_sum = prover.inner_prod(private_array.data(), public_array.data(), private_array.size());
	auto merkle_root_h = prover.commit_public_array(all_pub_mask, public_array.data(), log_length, inner_product_sum, all_sum);

	/* Verifier. */
	int proof_size;
	double v_time, p_time;
	auto processed = public_array_prepare_generic(public_array.data(), log_length);
	if (verifier.verify_poly_commitment(all_sum, log_length, processed, all_pub_mask, v_time, proof_size, p_time, merkle_root_l, merkle_root_h)) {
		cout << "Verification pass." << endl;
	}

    return 0;
}
