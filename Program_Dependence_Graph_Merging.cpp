#include <iostream>
#include <fstream>
#include <stdio.h>
#include <cstring>
#include <string>
#include <queue>
#include <algorithm>
#include <vector>
#include <stack>
#include <thread>
#include <mutex>
#include <stdlib.h>
#include <sstream>
#include <map>
#include <cmath> 
#include <time.h>

using namespace std;

class Edge {
public:
	int father_id; // from which graph
	string label;
	string style;
	int from_id, to_id; // the id of the Nodes
	Edge() {}
	Edge(int _father_id, string _label, string _style, int _from_id, int _to_id) {
		father_id = _father_id; label = _label; style = _style;
		from_id = _from_id; to_id = _to_id;
	}

	bool operator == (Edge &e) {
		if ((father_id == e.father_id) && (label == e.label) && (style == e.style)
			&& (from_id == e.from_id) && (to_id == e.to_id))
			return true;
		return false;
	}
};

class Node {
public:
	int id;
	int father_id; // from which graph
	int coord; // which line
	int al;			//whether it has been merged
	string label;
	string type, varname, vartype; // a groop of variates
	vector<Edge> itself_edges; //a circle ok?
	vector<Edge> out_edges; // edges start from the Node
	vector<Edge> in_edges; // edges ends at the node

	Node() {}
	Node(int _id, int _father_id, int _coord, string _label, string _type, string _varname, string _vartype) {
		id = _id; father_id = _father_id; coord = _coord;
		label = _label;
		type = _type; varname = _varname; vartype = _vartype;
		itself_edges.clear();
		in_edges.clear();
		out_edges.clear();
		al = 0;
	}

	bool operator == (Node &n1) {
		if (id == n1.id) {
			if (father_id == n1.father_id) {
				if (label == n1.label)
					return true;
			}
		}
		return false;
	}
};

class graph {
public:
	vector<Node> ENTRY_nodes;
	vector<Node> DECL_nodes;
	vector<Node> ASSIGN_nodes;
	vector<Node> CONTROL_nodes;
	vector<Node> CALL_nodes;
	vector<Node> RETURN_nodes;

	vector<Node> nodes; // store all the nodes, easy to code
	vector<Edge> edges;

	graph() {
		ENTRY_nodes.clear(); DECL_nodes.clear(); ASSIGN_nodes.clear();
		CONTROL_nodes.clear(); CALL_nodes.clear(); RETURN_nodes.clear();
		nodes.clear(); edges.clear();
	}
};

class merged_node {
public:
	vector<Node> nodes;
	int id;
	string type;

	merged_node() {
		nodes.clear();
		id = 0;
		type = "";
	}

	bool operator == (merged_node &mn) {
		if (type != mn.type) return false;
		if (nodes.size() != mn.nodes.size()) return false;
		vector<Node> ::iterator it1, it2;
		for (it1 = nodes.begin(); it1 != nodes.end(); ++it1) {
			int flag = 0;
			for (it2 = mn.nodes.begin(); it2 != mn.nodes.end(); ++it2) {
				if (*it1 == *it2) {
					flag = 1;
					break;
				}
			}
			if (flag == 0) return false;
		}
		return true;
	}
};

class merged_edge {
public:
	vector<Edge> edges;
	int from_id, to_id;
	string label;

	merged_edge() {
		edges.clear();
		from_id = 0; to_id = 0;
		label = "";
	}

	bool operator == (merged_edge &mm) {
		if (edges.size() != mm.edges.size()) return false;

		vector<Edge> ::iterator it1, it2;
		for (it1 = edges.begin(); it1 != edges.end(); ++it1) {
			int flag = 0;
			for (it2 = mm.edges.begin(); it2 != mm.edges.end(); ++it2) {
				if (*it1 == *it2) {
					flag = 1;
					break;
				}
			}
			if (flag == 0) return false;
		}
		return true;
	}
};

class individual {
public:
	vector<merged_node> m_nodes;
	vector<merged_edge> m_edges;
	individual() {
		m_nodes.clear();
		m_edges.clear();
	}
};

vector<graph> first_g;

graph load_graph(string FileName, int _father_id) {
	graph p;
	ifstream infile;
	infile.open(FileName.data());
	string s;

	vector<Node> _nodes;
	int hash_table[500];
	memset(hash_table, 0, sizeof(hash_table));
	int hash_ptr = 1;

	while (getline(infile, s)) {
		int temp_num = 0;
		string _label = "";
		_nodes.clear();

		//Node
		int _coord = 0;
		string _type = ""; string _varname = ""; string _vartype = "";

		//Edge
		string _style = "";
		int _from_id = 0, _to_id = 0;

		int ptr = 0; // the pointer of s
		if (s[0] != 'N') continue; //first line, contain no msg		
		for (ptr = 4; ; ++ptr) {
			if (s[ptr] >= '0' && s[ptr] <= '9') {
				temp_num *= 10;
				temp_num += s[ptr] - '0';
			}
			else break;
		}
		hash_table[hash_ptr] = temp_num;
		temp_num = hash_ptr;
		++hash_ptr;

		// get label
		ptr = s.find("label");
		if (ptr != s.npos) {
			ptr += 6;
			int l1 = s.find("\"", ptr);
			int l2 = s.find("\"", l1 + 1);
			_label = s.substr(l1 + 1, l2 - l1 - 1);
		}


		// judge if it is an edge
		ptr = s.find("->");
		if (ptr == s.npos) { // we don't find the "->", so it is a msg about Node
			ptr = s.find("type");
			if (ptr != s.npos) {
				ptr += 5;
				int l1 = s.find("\"", ptr);
				int l2 = s.find("\"", l1 + 1);
				_type = s.substr(l1 + 1, l2 - l1 - 1);
				if (_type == "CALL") {
					ptr = s.find("label");
					if (ptr != s.npos) {
						ptr += 6;
						int l1 = s.find("\"", ptr);
						int l2 = s.find(")", l1 + 1);
						_label = s.substr(l1 + 1, l2 - l1);
					}
				}
			}

			ptr = s.find("coord");
			if (ptr != s.npos) {
				ptr += 7;
				while (s[ptr] >= '0' && s[ptr] <= '9') {
					_coord *= 10;
					_coord += s[ptr] - '0';
					++ptr;
				}
			}

			ptr = s.find("varname");
			if (ptr != s.npos) {
				ptr += 9;
				int l1 = s.find("\"", ptr);
				int l2 = s.find("\"", l1 + 1);
				_varname = s.substr(l1 + 1, l2 - l1 - 1);
			}

			ptr = s.find("vartype");
			if (ptr != s.npos) {
				ptr += 9;
				int l1 = s.find("\"", ptr);
				int l2 = s.find("\"", l1 + 1);
				_vartype = s.substr(l1 + 1, l2 - l1 - 1);
			}

			Node temp_N(temp_num, _father_id, _coord, _label, _type, _varname, _vartype);

			p.nodes.push_back(temp_N);

		}
		else { // so it is a msg about Edge 
			for (ptr = 4; ; ++ptr) {
				if (s[ptr] >= '0' && s[ptr] <= '9') {
					_from_id *= 10;
					_from_id += s[ptr] - '0';
				}
				else break;
			}
			for (int i = 1; i < hash_ptr; ++i) {
				if (_from_id == hash_table[i]) {
					_from_id = i;
					break;
				}
			}

			ptr = s.find("->");
			ptr += 6;
			while (s[ptr] >= '0' && s[ptr] <= '9') {
				_to_id *= 10;
				_to_id += s[ptr] - '0';
				++ptr;
			}
			for (int i = 1; i < hash_ptr; ++i) {
				if (_to_id == hash_table[i]) {
					_to_id = i;
					break;
				}
			}

			ptr = s.find("style");
			if (ptr != s.npos) {
				if (s.find("dotted") != s.npos) _style = "dotted";
				else _style = "dashed";
			}


			Edge temp_edge(_father_id, _label, _style, _from_id, _to_id);
			p.edges.push_back(temp_edge);

			vector<Node> ::iterator it;
			for (it = p.nodes.begin(); it != p.nodes.end(); ++it) {
				if ((*it).id == temp_edge.from_id && (*it).id == temp_edge.to_id) {
					(*it).itself_edges.push_back(temp_edge);
					(*it).out_edges.push_back(temp_edge);
				}
				else if ((*it).id == temp_edge.from_id) {
					(*it).out_edges.push_back(temp_edge);
				}
				else if ((*it).id == temp_edge.to_id) { // use else if to avoid a edge appear twice
					(*it).in_edges.push_back(temp_edge);
				}
			}
		}
	}
	vector<Node> ::iterator it;
	for (it = p.nodes.begin(); it != p.nodes.end(); ++it) {
		if ((*it).type == "ENTRY") p.ENTRY_nodes.push_back(*it);
		else if ((*it).type == "DECL") p.DECL_nodes.push_back(*it);
		else if ((*it).type == "ASSIGN") p.ASSIGN_nodes.push_back(*it);
		else if ((*it).type == "CONTROL") p.CONTROL_nodes.push_back(*it);
		else if ((*it).type == "CALL") p.CALL_nodes.push_back(*it);
		else if ((*it).type == "RETURN") p.RETURN_nodes.push_back(*it);
	}
	return p;
}

// now solve the problem of output
void print_graph(individual ind, int filename) { // 0 = fibs, 1 = drug, 2 = odd
	ofstream fout;
	switch (filename) {
	case 0: fout.open("print_fibs.dot"); break;
	case 1: fout.open("print_drugs.dot"); break;
	case 2: fout.open("print_odd.dot"); break;
	}
	fout << "digraph PDG_Graph{" << endl;

	cout << "Nodes:" << ind.m_nodes.size() << endl;
	cout << "Edges:" << ind.m_edges.size() << endl;

	//now we need we delete the edges which are the same
	vector<merged_node> ::iterator it;
	vector<Node> ::iterator inner_it;

	for (it = ind.m_nodes.begin(); it != ind.m_nodes.end(); ++it) {
		inner_it = (*it).nodes.begin();

		if ((*it).type == "ENTRY") {
			fout << "Node" << (*it).id << " [type=\"ENTRY\",label=\"Entry\",shape=\"hexagon\"];" << endl;
			//only has one situation, so I make it a special occasion
		}

		else if ((*it).type == "DECL") {
			fout << "Node" << (*it).id << " [type=\"DECL\",label=\"" << (*it).id;
			for (inner_it = (*it).nodes.begin(); inner_it != (*it).nodes.end(); ++inner_it) {
				fout << "|" << (*inner_it).label;
			}
			fout << "\",shape=\"record\"];" << endl;
		}

		else if ((*it).type == "ASSIGN") {
			fout << "Node" << (*it).id << " [type=\"ASSIGN\",label=\"" << (*it).id;
			for (inner_it = (*it).nodes.begin(); inner_it != (*it).nodes.end(); ++inner_it) {
				fout << "|" << (*inner_it).label;
			}
			fout << "\",shape=\"Mrecord\"];" << endl;
		}

		else if ((*it).type == "CONTROL") {
			fout << "Node" << (*it).id << " [type=\"CONTROL\",label=\"" << (*it).id;
			for (inner_it = (*it).nodes.begin(); inner_it != (*it).nodes.end(); ++inner_it) {
				fout << "|" << (*inner_it).label;
			}
			fout << "\",shape=\"Mrecord\"];" << endl;
		}

		else if ((*it).type == "CALL") {
			fout << "Node" << (*it).id << " [type=\"CALL\",label=\"" << (*it).id;
			for (inner_it = (*it).nodes.begin(); inner_it != (*it).nodes.end(); ++inner_it) {
				fout << "|" << (*inner_it).label;
			}
			fout << "\",shape=\"Mrecord\"];" << endl;
		}

		else if ((*it).type == "RETURN") {
			fout << "Node" << (*it).id << " [type=\"RETURN\",label=\"" << (*it).id << "|return 0\"];" << endl;
		} // only one return
	}

	vector<merged_edge> ::iterator edge_it;
	for (edge_it = ind.m_edges.begin(); edge_it != ind.m_edges.end(); ++edge_it) {
		string _label = (*edge_it).label;
		fout << "Node" << (*edge_it).from_id << "->Node" << (*edge_it).to_id << " [label=\""
			<< _label << "\"";
		if (_label == "D" || _label == "LD") fout << ",style=dashed]" << endl;
		else if (_label == "DE" || _label == "DATA") fout << ",style=dotted]" << endl;
		else fout << "]" << endl;
	}
	fout << "}" << endl;
	fout.close(); // output ends
	cout << "RESULT GET" << endl;
}

void GetInfix(vector<string>& infix, string str) {
	infix.clear();
	string line;
	line = str;

	istringstream sin(line);
	string tmp;
	while (sin >> tmp)
	{
		infix.push_back(tmp);
	}
}

// 初始化操作符
void InitOperators(map<string, int>& opers) {
	opers.clear();
	opers["("] = 100;
	opers[")"] = 900;
	opers["+"] = 100;
	opers["-"] = 100;
	opers["*"] = 200;
	opers["/"] = 200;
}

bool IsOperator(const string& op, const map<string, int>& opers) {
	auto cit = opers.find(op);
	if (cit != opers.end())
	{
		return true;‘，                 
	}
	else
	{
		return false;
	}
}

void InfixToPrefix(const vector<string>& infix, vector<string>& prefix, map<string, int>& opers) {
	prefix.clear();
	stack<string> stk; // 操作符栈
	for (int i = infix.size() - 1; i >= 0; --i) // 从右到左扫描
	{
		if (!IsOperator(infix[i], opers)) // 如果是操作数
		{
			prefix.push_back(infix[i]);
		}
		else // 如果是操作符
		{
			if (infix[i] == ")") // 如果是右括号，则直接入栈
			{
				stk.push(infix[i]);
			}
			else if (infix[i] == "(") // 如果是左括号
			{
				// 依次弹出栈中的操作符，直至遇到右括号
				while (!stk.empty())
				{
					if (stk.top() == ")")
					{
						stk.pop();
						break;
					}
					else
					{
						prefix.push_back(stk.top());
						stk.pop();
					}
				}
			}
			else // 如果是其他操作符
			{
				while (!stk.empty() && stk.top() != ")" && opers[stk.top()] > opers[infix[i]]) // 栈顶操作符优先级大于当前操作符优先级
				{
					prefix.push_back(stk.top());
					stk.pop();
				}
				// 将当前操作符入栈
				stk.push(infix[i]);
			}
		}
	}

	// 检测操作符栈是否为空
	while (!stk.empty())
	{
		prefix.push_back(stk.top());
		stk.pop();
	}
	// 将prefix翻转
	reverse(prefix.begin(), prefix.end());
	return;
}

// get_pre 得到指定运算式的前缀表达式 
vector<string> get_pre(string ori) {
	map<string, int> opers;
	InitOperators(opers);

	vector<string> infix, prefix;
	GetInfix(infix, ori);
	InfixToPrefix(infix, prefix, opers);

	return prefix;
}

bool is_num(string str) {
	if (str[0] >= '0' && str[0] <= '9') return true;
	return false;
}

bool is_op(string str) {
	if (str == "+" || str == "-" || str == "*" || str == "/") return true;
	return false;
}

bool is_judge(string str) {
	if (str == ">" || str == ">=" || str == "<" || str == "<=" || str == "==") return true;
	return false;
}

int str2num(string str) {
	int n = 0, p = str.size();
	for (int i = 0; i < p; ++i) {
		n *= 10;
		n += str[i] - '0';
	}
	return n;
}

int med2(string &x, int m, string &y, int n) {
	int backup[25][25];
	//	int **backup = new int*[m];					//最后无法delete 可能成为问题？
	//	for (int i = 0; i < n; i++)				
	//		backup[i] = new int[n];
	memset(backup, 0, sizeof(backup));
	for (int i = 0; i <= n; ++i)
		backup[0][i] = i;
	for (int i = 0; i <= m; ++i)
		backup[i][0] = i;
	for (int i = 1; i <= m; ++i) {
		for (int j = 1; j <= n; ++j) {
			if (x[i - 1] == y[j - 1])
				backup[i][j] = backup[i - 1][j - 1];
			else {
				int a = backup[i - 1][j];
				int b = backup[i][j - 1];
				int c = backup[i - 1][j - 1];
				int temp;
				backup[i][j] = ((temp = a < b ? a : b) < c ? temp : c) + 1;
			}
		}
	}
	return backup[m][n];
}

double ASSIGN_AST(string x, string y) {
	double relation = 0;
	double struction = 0;

	int ptr_x = x.find("=");
	string ori_x = x.substr(0, ptr_x - 1);
	x = x.substr(ptr_x + 1, x.length());

	int ptr_y = y.find("=");
	string ori_y = y.substr(0, ptr_y - 1);
	y = y.substr(ptr_y + 1, y.length());

	vector<string> x_pre, y_pre;
	x_pre = get_pre(x);
	y_pre = get_pre(y);

	bool related_x = false, related_y = false;
	vector<string> ::iterator it;
	for (it = x_pre.begin(); it != x_pre.end(); ++it) {
		if (*it == ori_x) {
			related_x = true; break;
		}
	}
	for (it = y_pre.begin(); it != y_pre.end(); ++it) {
		if (*it == ori_y) {
			related_y = true; break;
		}
	}
	if (related_x == related_y) relation = 0.5;

	if (x_pre.size() == 1 && y_pre.size() == 1) {
		it = x_pre.begin();
		string tx = *it;
		it = y_pre.begin();
		string ty = *it;

		bool num_x = false, num_y = false;
		if (tx[0] >= '0' && tx[0] <= '9') num_x = true;
		if (ty[0] >= '0' && ty[0] <= '9') num_y = true;

		if ((!num_x) && (!num_y)) struction = 0.5;
		else if (num_x && num_y) {

			int k = 0;
			double sx = 0, sy = 0;
			while (k < tx.length()) {
				sx *= 10;
				sx += tx[k] - '0';
				++k;
			}
			k = 0;
			while (k < ty.length()) {
				sy *= 10;
				sy += ty[k] - '0';
				++k;
			}

			double result = sx / sy;
			if (result > 1) result = 1 / result;
			struction = result;
		}
		else struction = 0;

		return relation + struction;
	}

	string trans_x = "", trans_y = "";
	int x_size = x_pre.size(), y_size = y_pre.size();
	int len = x_size > y_size ? x_size : y_size;
	vector<string> ::iterator s_it;
	for (s_it = x_pre.begin(); s_it != x_pre.end(); ++s_it) {
		if (is_num(*s_it)) trans_x += "a";
		else if (is_op(*s_it)) trans_x += "b";
		else trans_x += "c";
	}
	for (s_it = y_pre.begin(); s_it != y_pre.end(); ++s_it) {
		if (is_num(*s_it)) trans_y += "a";
		else if (is_op(*s_it)) trans_y += "b";
		else trans_y += "c";
	}

	int res = med2(trans_x, x_size, trans_y, y_size);
	struction = 0.5 * ((len - res) / len);

	return relation + struction;
	// whether related to itself take a weight of 50%
	// the other 50% is decided by the sentence 
}

double single_CONTROL(string x, string y) { // single control word's similarity
	double w_op = 1, w_num = 1;

	bool x_l = ((x.find("<") != x.npos) && (x.find("<=") == x.npos));
	bool x_le = (x.find("<=") != x.npos);
	bool x_g = ((x.find(">") != x.npos) && (x.find(">=") != x.npos));
	bool x_ge = (x.find(">=") != x.npos);
	bool x_e = ((x.find("==") != x.npos) || (x.find("!=") != x.npos) || ((x.find("=") != x.npos) && (!x_le) && (!x_ge)));

	bool y_l = ((y.find("<") != y.npos) && (y.find("<=") == y.npos));
	bool y_le = (y.find("<=") != y.npos);
	bool y_g = ((y.find(">") != y.npos) && (y.find(">=") != y.npos));
	bool y_ge = (y.find(">=") != y.npos);
	bool y_e = ((y.find("==") != y.npos) || (y.find("!=") != y.npos) || ((y.find("=") != y.npos) && (!y_le) && (!y_ge)));

	vector<string> x_pre = get_pre(x), y_pre = get_pre(y);
	vector<string> ::iterator x_it, y_it;

	if (x_e && y_e) {
		bool x_num = false, y_num = false;
		for (x_it = x_pre.begin(); x_it != x_pre.end(); ++x_it) {
			if (is_num(*x_it)) {
				x_num = true; break;
			}
		}
		for (y_it = y_pre.begin(); y_it != y_pre.end(); ++y_it) {
			if (is_num(*y_it)) {
				y_num = true; break;
			}
		}
		if (x_num == y_num) return 1;
		else return 0;
	}

	if ((x_e && (!y_e)) || ((!x_e) && y_e)) {
		return 0;
	}

	if ((x_l == y_l) && (x_le == y_le) && (x_g == y_g) && (x_ge == y_ge))
		w_op = 1;
	else if ((x_l == y_le) && (x_le == y_l) && (x_g == y_g) && (x_ge == y_ge))
		w_op = 0.9;
	else if ((x_l == y_l) && (x_le == y_le) && (x_ge == y_g) && (x_g == y_ge))
		w_op = 0.9;
	else if ((x_l == y_g) && (x_le == y_le) && (x_g == y_l) && (x_ge == y_ge))
		w_op = 0.4;
	else if ((x_l == y_l) && (x_le == y_ge) && (x_g == y_g) && (x_ge == y_le))
		w_op = 0.4;
	else if ((x_l == y_ge) && (x_le == y_le) && (x_g == y_g) && (x_ge == y_l))
		w_op = 0.3;
	else if ((x_l == y_l) && (x_g == y_le) && (x_le == y_g) && (x_ge == y_ge))
		w_op = 0.3;

	bool x_num = false, y_num = false;
	bool x_flag = false, y_flag = false;
	double xx = 0, yy = 0;

	for (x_it = x_pre.begin(); x_it != x_pre.end(); ++x_it) {
		if (is_num(*x_it)) {
			xx = str2num(*x_it);
			x_flag = true;
			break;
		}
	}

	for (y_it = y_pre.begin(); y_it != y_pre.end(); ++y_it) {
		if (is_num(*y_it)) {
			yy = str2num(*y_it);
			y_flag = true;
			break;
		}
	}

	if ((!x_flag) && (!y_flag)) return w_op;
	else if (((!x_flag) && y_flag) || x_flag && (!y_flag)) return w_op * 0.2;
	else {
		w_num = xx / yy;
		if (w_num > 1) w_num = yy / xx;
		return w_op * w_num;
	}
}

double CONTROL_AST(string x, string y) {
	double w_and = 1;

	bool and_x = false, and_y = false;
	and_x = (x.find("&&") != x.npos);
	and_y = (y.find("&&") != y.npos);

	//	if (and_x != and_y) w_and = 0.75;

	if (and_x && and_y) {
		int ptr_x = x.find("&&"), ptr_y = y.find("&&");
		string x1 = x.substr(0, ptr_x);
		string x2 = x.substr(ptr_x + 2, x.size());
		string y1 = y.substr(0, ptr_y);
		string y2 = y.substr(ptr_y + 2, y.size());

		double r1 = 0.5 * single_CONTROL(x1, y1) + 0.5 * single_CONTROL(x2, y2);
		double r2 = 0.5 * single_CONTROL(x1, y2) + 0.5 * single_CONTROL(x2, y1);
		double r = r1 > r2 ? r1 : r2;
		return w_and * r;
	}

	else if (and_x && (!and_y)) {
		w_and = 0.75;
		int ptr_x = x.find("&&");
		string x1 = x.substr(0, ptr_x);
		string x2 = x.substr(ptr_x + 2, x.size());

		double r1 = single_CONTROL(x1, y);
		double r2 = single_CONTROL(x2, y);
		double r = r1 > r2 ? r1 : r2;
		return w_and * r;
	}

	else if ((!and_x) && and_y) {
		w_and = 0.75;
		int ptr_y = y.find("&&");
		string y1 = y.substr(0, ptr_y);
		string y2 = y.substr(ptr_y + 2, y.size());

		double r1 = single_CONTROL(x, y1);
		double r2 = single_CONTROL(x, y2);
		double r = r1 > r2 ? r1 : r2;
		return w_and * r;
	}

	else if ((!and_x) && (!and_y)) {
		double r = single_CONTROL(x, y);
		return r;
	}

	return 0;
}

double Node_sim(Node x, Node y) {

	//first get all the special occasions solved
	if (x.type != y.type) return 0;

	if (x.type == "ENTRY") return 1.0;
	if (x.type == "RETURN") return 1.0;

	if (x.label == "Entry" || x.label == "ENTRY") {
		if (y.label == "Entry" || x.label == "ENTRY") return 1.0;
		else return 0;
	}
	if (x.label.find("break") != x.label.npos) {
		if (y.label.find("break") != y.label.npos) return 1.0;
		else return 0;
	}
	if (x.label.find("continue") != x.label.npos) {
		if (y.label.find("continue") != y.label.npos) return 1.0;
		else return 0;
	}
	if (x.label.find("return 0") != x.label.npos) {
		if (y.label.find("return 0") != y.label.npos) return 1.0;
		else return 0;
	}
	if (x.label.find("int main") != x.label.npos) {
		if (y.label.find("int main") != y.label.npos) return 1.0;
		else return 0;
	}
	// now all the special occasions are solved
	// can be added at any time

	// let's start the universial situations
	if (x.type == "DECL") {
		if (x.vartype == y.vartype)
			return 1;
		else if (x.vartype[0] != y.vartype[0])
			return 0;
		else if (x.vartype.find("*") != x.vartype.npos || x.vartype.find("[]") != x.vartype.npos)		//类型相同维度不同 
			if (y.vartype.find("*") != y.vartype.npos || y.vartype.find("[]") != y.vartype.npos)
				return 0.75;
		return 0.5;
	}

	if (x.type == "CALL") {
		if (x.label[0] == y.label[0]) { // are both input or output
			string x_var = "", y_var = "";
			string x_temp = x.label, y_temp = y.label;
			bool x_array = false, y_array = false;

			int l1 = x_temp.find("%", 0);
			int l2 = x_temp.find("\\", l1 + 1);
			x_var = x_temp.substr(l1 + 1, l2 - l1 - 1);
			if (x_temp.find("[") != x_temp.npos) x_array = true;

			l1 = y_temp.find("%", 0);
			l2 = y_temp.find("\\", l1 + 1);
			y_var = y_temp.substr(l1 + 1, l2 - l1 - 1);
			if (y_temp.find("[") != y_temp.npos) y_array = true;

			if (x_var == y_var && x_array == y_array) return 1;
			else if (x_var != y_var) return 0;
			else if (x_var == y_var && x_array != y_array) return 0.5;
		}
		else return 0;
	}

	if (x.type == "ASSIGN") {
		return ASSIGN_AST(x.label, y.label);
	}

	if (x.type == "CONTROL") {
		return CONTROL_AST(x.label, y.label);
	}

	return 0;
}

individual nodes2graph(vector<merged_node> v) {
	individual ind;
	ind.m_nodes = v;

	vector<merged_node> ::iterator it;
	vector<Node> ::iterator inner_it;
	vector<Edge> ::iterator edge_it;

	vector<merged_edge> ::iterator EDGE_it;

	int hash_table[15][100];
	memset(hash_table, 0, sizeof(hash_table));
	for (it = ind.m_nodes.begin(); it != ind.m_nodes.end(); ++it) {
		(*it).id = (it - ind.m_nodes.begin()) + 1;
		for (inner_it = (*it).nodes.begin(); inner_it != (*it).nodes.end(); ++inner_it) {
			hash_table[(*inner_it).father_id][(*inner_it).id] = (*it).id;
		}
	}

	for (it = ind.m_nodes.begin(); it != ind.m_nodes.end(); ++it) {
		for (inner_it = (*it).nodes.begin(); inner_it != (*it).nodes.end(); ++inner_it) {
			for (edge_it = (*inner_it).out_edges.begin(); edge_it != (*inner_it).out_edges.end(); ++edge_it) {
				int _from_id = hash_table[(*edge_it).father_id][(*edge_it).from_id];
				int _to_id = hash_table[(*edge_it).father_id][(*edge_it).to_id];
				string _label = (*edge_it).label;
				Edge temp = *edge_it;

				merged_edge TEMP;
				TEMP.edges.push_back(temp);
				TEMP.from_id = _from_id;
				TEMP.to_id = _to_id;
				TEMP.label = _label;

				bool flag = true;
				for (EDGE_it = ind.m_edges.begin(); EDGE_it != ind.m_edges.end(); ++EDGE_it) {
					if ((*EDGE_it).from_id == _from_id && (*EDGE_it).to_id == _to_id && (*EDGE_it).label == _label) {
						(*EDGE_it).edges.push_back(temp);
						flag = false;
						break;
					}
				}
				if (flag) {
					ind.m_edges.push_back(TEMP);
				}
			}
		}
	}

	return ind;
}

const int Max = 9999;      //假设一个图最多有这么多节点
						   //要求移动路径最小，就是两个融合图各自保留的节点最多,转换为带权二部图匹配问题

int can(int t, int m, int* link, double* visx, double* visy, double *lx, double *ly, double weight[200][200])				//检查是否能够生成增广路径  
{
	visx[t] = 1;
	for (int i = 1; i <= m; i++) {
		if (!visy[i] && lx[t] + ly[i] == weight[t][i]) {					//这里“lx[t]+ly[i]==w[t][i]”决定了这是在相等子图中找增广路的前提，非常重要
			visy[i] = 1;
			if (link[i] == -1 || can(link[i], m, link, visx, visy, lx, ly, weight))
			{
				link[i] = t;
				return 1;
			}
		}
	}
	return 0;
}

vector<merged_node> MaximumMatching(vector<merged_node> x, vector<Node> y)
{
	vector<Node> yy;
	yy = y;
	//			cout << "graph=" << y[0].father_id<<endl;
	//			cout << "x.size" << x.size() << endl;
	//			cout << "y.size" << y.size() << endl;
	double weight[200][200];		   //匹配权重
	memset(weight, 0, sizeof(weight));
	double visx[200], visy[200];
	double lx[200], ly[200];      //lx为x中的标杆值，ly为y中的标杆值
	int n = x.size();
	int m = y.size();
	int link[499];		   //link代表当前与Y集合中配对的X集合中的点
	memset(link, -1, sizeof(link));
	for (int i = 0; i < n; i++)				//计算匹配权重
		for (int j = 0; j < m; j++)
		{
			for (int ii = 0; ii<x[i].nodes.size(); ii++)
			{
				weight[i][j] += Node_sim(x[i].nodes[ii], y[j]);
			}
		}
	memset(ly, 0, sizeof(ly));				//初始y的标杆值为0
	for (int i = 0; i < n; i++)			//初始x的标杆值为x所能匹配的最大值
	{
		lx[i] = -1;
		for (int j = 0; j < m; j++)
		{
			if (lx[i] < weight[i][j])
				lx[i] = weight[i][j];
		}
	}
	for (int i = 0; i < n; i++)
	{
		int maxnum = 100;
		while (maxnum>0)
		{
			maxnum--;
			//			cout << "maxnum=" << maxnum << endl;
			memset(visx, 0, sizeof(visx));
			memset(visy, 0, sizeof(visy));
			if (can(i, m, link, visx, visy, lx, ly, weight))//如果它能够形成一条增广路径，那么就break
				break;
			//找新匹配点时，能够进入二分子图的边的条件为L(x) + L(y) >= weight(xy)
			double d = -1;//后面应该加入新的边,这里应该先计算d值
			for (int j = 0; j < n; j++)//对于搜索过的路径上的XY点，设该路径上的X顶点集为S，Y顶点集为T，对所有在S中的点xi及不在T中的点yj
				if (visx[j])
					for (int k = 0; k < m; k++)
						if (!visy[k])
						{
							if (d < lx[j] + ly[k] - weight[j][k])
								d = lx[j] + ly[k] - weight[j][k];
						}
			if (d == -1)					//找不到可以加入的边，找不到完美匹配
				break;
			for (int j = 0; j < n; j++)	//更新标杆值
				if (visx[j])
					lx[j] -= d;
			for (int j = 0; j < m; j++)
				if (visy[j])
					ly[j] += d;
		}
	}
	int n_n = n + 1;
	int flage = 1;
	int flagr = 1;
	int flagi = 1;
	for (int i = 0; i < n; i++)
	{
		if (flage)
		{
			if (x[i].type == "ENTRY")
			{
				for (int j = 0; j < m; j++)
					if (yy[j].type == "ENTRY")
					{
						x[i].nodes.push_back(yy[j]);
						flage--;
						link[i] = -1;
						break;
					}
				link[i] = -1;
				i++;
			}
		}
		if (i >= n)
			break;
		if (flagr)
		{
			if (x[i].nodes[0].label == "RETURN 0")
			{
				for (int j = m - 1; j >= 0; j--)
				{
					if (yy[j].label == "RETURN 0")
					{
						x[i].nodes.push_back(yy[j]);
						flagr--;
						break;
					}
				}
				link[i] = -1;
				i++;
			}
		}
		if (i >= n)
			break;
		if (flagi)
		{
			if (x[i].nodes[0].label == "int main()")
			{
				for (int j = m - 1; j >= 0; j--)
					if (yy[j].label == "int main()")
					{
						x[i].nodes.push_back(yy[j]);
						flagi--;
						link[i] = -1;
						break;
					}
				link[i] = -1;
				i++;
			}
		}
		if (i >= n)
			break;
		if (link[i] != -1 && link[i]<m&&link[i] >= 0)
		{
			if (x[i].type != yy[link[i]].type)
				continue;
			x[i].nodes.push_back(yy[link[i]]);
			link[i] = -1;
		}
	}
	int time = 3;
	for (int i = 0; i < m; i++)
		if (link[i] != -1 && link[i] < m&&link[i] >= 0)
		{
			if (time)
			{
				time--;
				continue;
			}
			time = 3;
			merged_node temp;
			temp.nodes.push_back(yy[link[i]]);
			//			cout << "bug1";
			temp.id = n_n;
			//			cout << "bug2";
			temp.type = yy[link[i]].type;
			//			cout << "bug3";
			x.push_back(temp);
			n_n++;
			//			cout << "not here";
		}
	int mmax = x.size();
	int flag1 = 0;
	int flag11 = 0;
	for (int i = 0; i < mmax; i++)
	{
		if (x[i].nodes[0].label == "return 0")
		{
			if (flag1 == 0)
			{
				flag1++;
				flag11 = i;
			}
			else
			{
				x[i].nodes.insert(x[i].nodes.end(), x[flag11].nodes.begin(), x[flag11].nodes.end());
				x.erase(x.begin() + flag11);
				mmax--;
			}

		}
	}
	mmax = x.size();
	flag1 = 0;
	flag11 = 0;
	for (int i = 0; i < mmax; i++)
	{
		if (x[i].nodes[0].label == "int main()")
		{
			if (flag1 == 0)
			{
				flag1++;
				flag11 = i;
			}
			else
			{
				x[i].nodes.insert(x[i].nodes.end(), x[flag11].nodes.begin(), x[flag11].nodes.end());
				x.erase(x.begin()+flag11);
				mmax--;
			}

		}
	}
	flag1 = 0;
	flag11 = 0;
	for (int i = 0; i < mmax; i++)
	{
		if (x[i].type == "ENTRY")
		{
			if (flag1 == 0)
			{
				flag1++;
				flag11 = i;
			}
			else
			{
				x[i].nodes.insert(x[i].nodes.end(), x[flag11].nodes.begin(), x[flag11].nodes.end());
				x.erase(x.begin() + flag11);
				mmax--;
			}

		}
	}
	return x;
}

int main() {
	srand((unsigned)time(NULL));
	int input;
	cin >> input;

	if (input == 0) {
		first_g.push_back(load_graph("FIBS01.dot", 0));
		first_g.push_back(load_graph("FIBS02.dot", 1));
		first_g.push_back(load_graph("FIBS03.dot", 2));
	}
	else if (input == 1) {
		first_g.push_back(load_graph("DRUG01.dot", 0));
		first_g.push_back(load_graph("DRUG02.dot", 1));
		first_g.push_back(load_graph("DRUG03.dot", 2));
		first_g.push_back(load_graph("DRUG04.dot", 3));
		first_g.push_back(load_graph("DRUG05.dot", 4));
		first_g.push_back(load_graph("DRUG06.dot", 5));
		first_g.push_back(load_graph("DRUG07.dot", 6));
		first_g.push_back(load_graph("DRUG08.dot", 7));
		first_g.push_back(load_graph("DRUG09.dot", 8));
		first_g.push_back(load_graph("DRUG10.dot", 9));
	}
	else if (input == 2) {
		first_g.push_back(load_graph("ODD01.dot", 0));
		first_g.push_back(load_graph("ODD02.dot", 1));
		first_g.push_back(load_graph("ODD03.dot", 2));
		first_g.push_back(load_graph("ODD04.dot", 3));
		first_g.push_back(load_graph("ODD05.dot", 4));
		first_g.push_back(load_graph("ODD06.dot", 5));
		first_g.push_back(load_graph("ODD07.dot", 6));
		first_g.push_back(load_graph("ODD08.dot", 7));
		first_g.push_back(load_graph("ODD09.dot", 8));
		first_g.push_back(load_graph("ODD10.dot", 9));
	}

	if (input == 0) {
		individual ind, ff;
		for (int ii = 0; ii < 50; ++ii) { //随机顺序融合以求得到更优解 
			ind.m_nodes.clear();
			ind.m_edges.clear();
			int temp[10], re[10];
			memset(temp, 0, sizeof(temp));
			memset(re, 0, sizeof(re));
			for (int j = 0; j < 3; ++j) {
				while (1) {
					int t = rand() % 3;
					if (temp[t] == 0) {
						temp[t] = 1;
						re[j] = t;
						break;
					}
				}
			}
			int fst = re[0];
			vector<merged_node> temp_m;
			vector<Node> ::iterator it;
			graph g_temp;
			g_temp = first_g[fst];
			for (it = g_temp.nodes.begin(); it != g_temp.nodes.end(); ++it) {
				merged_node r;
				r.id = it - g_temp.nodes.begin() + 1;
				r.nodes.push_back(*it);
				r.type = (*it).type;
				ind.m_nodes.push_back(r);
			}
			for (int i = 1; i < 3; ++i) {
				vector<merged_node> t1 = ind.m_nodes;
				vector<Node> t2 = (first_g[re[i]]).nodes;
				ind.m_nodes = MaximumMatching(t1, t2);
			}
			if (ii != 0) {
				if (ff.m_edges.size() >= ind.m_edges.size() && ff.m_nodes.size() >= ind.m_nodes.size()) {
					ff = ind;
				}
			}
			else {
				ff = ind;
			}
			//compare 
		}
		ff = nodes2graph(ff.m_nodes);
		print_graph(ff, input);
	}
	else {
		individual ind, ff;
		for (int ii = 0; ii < 50; ++ii) { //随机顺序融合以求得到更优解 
			ind.m_nodes.clear();
			ind.m_edges.clear();
			int temp[15], re[15];
			memset(temp, 0, sizeof(temp));
			memset(re, 0, sizeof(re));
			for (int j = 0; j < 10; ++j) {
				while (1) {
					int t = rand() % 10;
					if (temp[t] == 0) {
						temp[t] = 1;
						re[j] = t;
						break;
					}
				}
			}

			int fst = re[0];
			vector<merged_node> temp_m;
			vector<Node> ::iterator it;
			graph g_temp;
			g_temp = first_g[fst];
			for (it = g_temp.nodes.begin(); it != g_temp.nodes.end(); ++it) {
				merged_node r;
				r.id = it - g_temp.nodes.begin() + 1;
				r.nodes.push_back(*it);
				r.type = (*it).type;
				ind.m_nodes.push_back(r);
			}
			for (int i = 1; i < 10; ++i) {
				vector<merged_node> t1 = ind.m_nodes;
				vector<Node> t2 = (first_g[re[i]]).nodes;
				ind.m_nodes = MaximumMatching(t1, t2);
			}
			cout << endl;
			if (ii != 0) {
				if (ff.m_edges.size() >= ind.m_edges.size() && ff.m_nodes.size() >= ind.m_nodes.size()) {
					ff = ind;
				}
			}
			else {
				ff = ind;
			}
			//compare 
		}
		ff = nodes2graph(ff.m_nodes);
		print_graph(ff, input);
	}
	system("pause");
	return 0;
}



