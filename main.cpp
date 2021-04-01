#include <iostream>
#include <vector>
#include <cmath>
#include <map>
#include <ctime>
#include "BigInt.hpp"
#include <fstream>

#define VAR 3
using namespace std;
long s_time;
time_t begin_time, end_time;
int try_count = 10;
string word = "Well I'm not here cause I wanna be here Just dunno anything else";/*sSpent so long deep in love with you Fell out of love with myself I'm screaming at your memory I'm saying  get the hell right out of my house So I take your pictures off the wall And it's just me and Johnny Walker left around You can lie lie lie You can lie to me I don't care I'm not losing sleep If it don't bother you It sure as hell don't bother meome words I never say to you my darling */
vector<long> simples;

void readSimple() {
    fstream in("simple.txt");
    if(in.is_open()) {
        in.seekg(ios_base::beg);
        long s;
        for (in >> s; !in.eof(); in >> s) {
            simples.push_back(s);
        }
    }
    in.close();
}

void writeSimples(long p, long g) {
    fstream out("simple.txt", ios::out|ios::in|ios::trunc);
    simples.push_back(p);
    simples.push_back(g);
    if(out.is_open()) {
        for (auto artist : simples) {
            out << artist << endl;
        }
    }
    out.close();
}

void writeTimes(double t, const BigInt& n, const BigInt& e, const BigInt& d) {
    fstream out("times.txt", ios::out|ios::app);

    if(out.is_open()) {
#if VAR == 1
        out << "TD ";
#elif VAR == 2
        out << "MR ";
#elif VAR == 3
        out << "SS ";
#endif
            out << t  << ' ' << try_count << ' ' << word.size()  << ' ' << n<< ' ' << e<< ' ' << d << endl;
    }
    out.close();
}

int Euclid(int A, int B) {
    while (A % B != 0) {
        int C = A % B;
        A = B;
        B = C;
    }
    return B;
}

BigInt Euclid(BigInt A, BigInt B) {
    while (A % B != 0) {
        BigInt C = A % B;
        A = B;
        B = C;
    }
    return B;
}

BigInt extendedEuclid(BigInt A, BigInt B) {
    vector<BigInt> mod;
    vector<BigInt> div;
    int size = 0;
    do {
        mod.push_back(A % B);
        div.push_back((BigInt) A / B);
        A = B;
        B = mod[size];
        size++;
    } while (A % B != 0);
    BigInt x = 0, y = 0, x1 = 0, y1 = 1;
    for (int i = size - 1; i >= 0; --i) {
        x = y1;
        y = x1 - y1 * div[i];
        x1 = x;
        y1 = y;
    }
    return y;
}

long extendedEuclid(long A, long B) {
    vector<long> mod;
    vector<long> div;
    int size = 0;
    do {
        mod.push_back(A % B);
        div.push_back((long) A / B);
        A = B;
        B = mod[size];
        size++;
    } while (A % B != 0);
    long x = 0, y = 0, x1 = 0, y1 = 1;
    for (int i = size - 1; i >= 0; --i) {
        x = y1;
        y = x1 - y1 * div[i];
        x1 = x;
        y1 = y;
    }
    return y;
}

long powm(long a, long b, long n) {
    long c = 1;
    while (b) {
        if(b % 2 == 0) {
            b /= 2;
            a = (a * a) % n;
        } else {
            b--;
            c = (c * a) % n;
        }
    }
    return c;
}

map<BigInt, BigInt> factorize(BigInt x) {
    map<BigInt, BigInt> factors;
    auto g = sqrt(x);
    for (int i = 2; i <= sqrt(x); i++) {
        while (x % i == 0) {
            factors[i]++;
            x /= i;
        }
    }
    if(x != 1) {
        factors[x]++;
    }
    return factors;
}

BigInt powm(BigInt a, BigInt b, const BigInt &n) {
    BigInt c = 1;
    while (b > 0) {
        if(b % 2 == 0) {
            b /= 2;
            a = (a * a) % n;
        } else {
            b--;
            c = (c * a) % n;
        }
    }
    return c;
}

int trialDivision(int n) {
    begin_time = time(nullptr);
    for (int i = 2; i < floor(sqrt(n)); i++) {
        if(n % i == 0) {
            return i;
        }
    }
    end_time = time(nullptr);
    cout << n << " is prime TD for " << end_time - begin_time << " sec" << endl;
    return 0;
}

BigInt trialDivision(const BigInt &n) {
    auto sqrtn = sqrt(n);
    begin_time = time(nullptr);
    BigInt max = 0;
    for (BigInt i = 2; i < sqrtn; i++) {
        if(i > max) {
            max = i;
        }else {
            exit(6);
        }
        if(n % i == 0) {
            return i;
        }
    }
    end_time = time(nullptr);
    return 0;
}

pair<long, BigInt> getSD(const BigInt &n) {
    long s = -1;
    BigInt d = 1;
    auto fact_seq = factorize(n);
    for (const auto &a: fact_seq) {
        if(a.first == 2) {
            s = a.second.to_long();
        } else {
            d *= pow(a.first, a.second);
        }
    }
    if(d % 2 == 0) {
        cout << "screwed up" << endl;
    }
    return make_pair(s, d);
}

bool stepWithSeq(const BigInt &x0, const BigInt &n, long s) {
    BigInt xi = x0, xi1;
    for (int i = 0; i < s - 1; ++i) {
        xi1 = powm(xi, 2, n);
        xi = xi1;
        if(n - 1 == xi1) {
            return true;
        }
    }
    return false;
}

BigInt getRandom(const BigInt &min, const BigInt &max, long random_num = 1) {
    auto buff = s_time;
    s_time = time(nullptr) + min.to_int() + random_num;
    if(s_time == buff)
        s_time -= time(nullptr) / (min.to_int() + 1);
    srand(s_time);
    BigInt random_variable = rand() % max + min;
    if(random_variable > max - 1) {
        random_variable = max - 2;
    }
    return random_variable;
}

bool primeMillerRabinTest(const BigInt &n, int r) {
    begin_time = time(nullptr);
    BigInt x0;
    if(n > 2 && n % 2 != 0) {
        int witness = 0;
        auto fact_seq = factorize(n - 1);
        auto s_d = getSD(n - 1);
        for (int a = 2; a < r + 1; ++a) {
            x0 = powm(a, s_d.second, n);
            if(x0 == 1 || x0 == n - 1) {
                witness++;
                continue;
            } else {
                if(stepWithSeq(x0, n, s_d.first)) {
                    witness++;
                    continue;
                } else {
                    break;
                }
            }

        }
        if(witness == r - 1) {
            end_time = time(nullptr);
            return true;
        }
    }
    return false;
}

////////////////////////////////////Legendre!!//////////////////////////////
void trowingOutEvenA(map<pair<BigInt, BigInt>, BigInt> &fact_seq) {
    vector<pair<BigInt, BigInt>> to_delete;
    for (const auto &a : fact_seq) {
        if(a.second % 2 == 0) {
            to_delete.push_back(a.first);
        }
    }
    for (const auto &e: to_delete) {
        fact_seq.erase(e);
    }
}

int findEasyOnes(map<pair<BigInt, BigInt>, BigInt> &fact_seq) {
    vector<pair<BigInt, BigInt>> to_delete;
    int res = 1;
    BigInt sign = 0;
    BigInt one = -1;
    for (const auto &a : fact_seq) {
        if(a.first.first == 2) {
            auto power = (pow(a.first.second, 2) - 1) / 8;
            sign = pow(one, power);
            res *= sign.to_int();
            to_delete.push_back(a.first);
        }
    }
    for (const auto &e: to_delete) {
        fact_seq.erase(e);
    }
    return res;
}

int getTrickyOnes(map<pair<BigInt, BigInt>, BigInt> &fact_seq) {
    vector<pair<BigInt, BigInt>> to_delete;
    vector<pair<pair<BigInt, BigInt>, BigInt>> to_insert;
    int res = 1;
    BigInt sign = 0;
    BigInt one = -1;
    for (const auto &a : fact_seq) {
        auto power = ((a.first.first - 1) / 2) * ((a.first.second - 1) / 2);
        sign = pow(one, power);
        res *= sign.to_int();
        to_insert.emplace_back(a);
        to_delete.push_back(a.first);
    }
    for (const auto &e: to_delete) {
        fact_seq.erase(e);
    }
    for (const auto &a: to_insert) {
        fact_seq.emplace(make_pair(a.first.second, a.first.first), a.second);
    }
    return res;
}

map<pair<BigInt, BigInt>, BigInt> decomposition(const map<BigInt, BigInt> &fact_seq, const BigInt &p) {
    map<pair<BigInt, BigInt>, BigInt> dec;
    for (const auto &ai: fact_seq) {
        dec.emplace(make_pair(ai.first, p), ai.second);
    }
    return dec;
}

int getNegativeOnes(BigInt &a, const BigInt &p) {
    int res = 1;
    BigInt sign = 0;
    if(a < 0) {
        sign = powm(-1, (p - 1) / 2, p);
        res *= sign.to_int();
        a = -a;
    }
    return res;
}

int LegendreSymbol(BigInt a, const BigInt &p) {
    int res = 1;
    auto neg_sign = getNegativeOnes(a, p);
    a = a % p;
    auto fact = factorize(a);
    auto dec = decomposition(fact, p);
    trowingOutEvenA(dec);
    auto sign = findEasyOnes(dec);
    auto sec_sign = getTrickyOnes(dec);
    res *= sec_sign * sign * neg_sign;
    if(dec.empty()) {
        return res;
    }
    for (const auto &e: dec) {
        res *= LegendreSymbol(e.first.first, e.first.second);
    }
    return res;
}
///////////////////////////////////Legendre!!//////////////////////////////

BigInt binpow(const BigInt &a, const BigInt &n) {
    if(n == 0)
        return 1;
    if(n % 2 == 1)
        return binpow(a, n - 1) * a;
    else {
        BigInt b = binpow(a, n / 2);
        return b * b;
    }
}

BigInt fastpow(BigInt value, BigInt pow) {
    BigInt result = 1;
    while (pow > 0) {
        if(pow % 2 == 1) {
            result *= value;
        }
        value *= value;
        pow /= 2;
    }
    return result;
}

bool primeSolovayStrassenTest(const BigInt &n, int k) {
    begin_time = time(nullptr);
    for (int i = 0; i < k; ++i) {
        BigInt a = getRandom(3, n);
        auto d = Euclid(a, n);
        if(d > 1) {
            return false;
        } else {
            auto left = LegendreSymbol(a, n);
            auto right = powm(a, (n - 1) / 2, n);
            if(left == -1) {
                if(right + 1 == n) {
                    continue;
                }
            }
            if(left != right) {
                return false;
            }
        }
    }
    end_time = time(nullptr);
    return true;
}

BigInt keyGeneration(int attempts, long random_num = 1, bool out = true) {
    auto total_time = time(nullptr);
    BigInt n;
    int i = 0;
    do {
        unsigned int num = getRandom(0, simples.size(), random_num + i).to_int();
        if(num > simples.size()) {
            exit(1);
        }
        auto p = simples[num];

        i++;
        auto R = getRandom(p, 4 * p + 2, random_num + i);
        if(R % 2 != 0) {
            R++;
        }
        n = p * R + 1;
#if VAR == 3
        }while(!primeSolovayStrassenTest(n, attempts));
    if(out)
        cout<< "SS: " ;
#elif VAR == 2
    } while (!primeMillerRabinTest(n, attempts));
    if(out)
        cout << "MR: ";
#elif VAR == 1
    } while (trialDivision(n) != 0);
    if(out)
        cout<< "TD: " ;
#endif
    if(out) {
        cout << n << " is prime for prime search " << difftime(time(nullptr), begin_time)/*end_time - begin_time*/
             << " sec" << endl;
        cout << "Total time with " << attempts << " checks is "
             << difftime(time(nullptr), total_time) /*time(nullptr) - total_time */<< '\n';
    }
    return n;
}

BigInt euler(const BigInt &p, const BigInt &q) {
    return (p - 1) * (q - 1);
}

BigInt enc(const BigInt &c, const BigInt &e, const BigInt &n) {
    return powm(c, e, n);
}

vector<BigInt> encryption(const string &some_word, const BigInt &e, const BigInt &n) {
    vector<BigInt> cypher;
    for (auto c : some_word) {
        auto h = enc(c, e, n);
        cypher.push_back(h);
    }
    return cypher;
}

BigInt dec(const BigInt &h, const BigInt &d, const BigInt &n) {
    return powm(h, d, n);
}

string decryption(const vector<BigInt> &code, const BigInt &d, const BigInt &n) {
    string decipher;
    for (const auto &h : code) {
        auto c = dec(h, d, n);
        decipher += c.to_int();
    }
    return decipher;
}

int main() {
    readSimple();
    auto start = time(nullptr);
    auto g = keyGeneration(try_count);
    auto halfg = g/2;
    if(LONG_MAX < halfg) {
        halfg = LONG_MAX/2;
    }
    auto p = keyGeneration(try_count, halfg.to_long());
    auto halfp = p/2;
    if(LONG_MAX < halfp) {
        halfp = LONG_MAX/2;
    }
    int i = 0;
    while (p == g) {
        i++;
        p = keyGeneration(try_count, i);
    }
    auto n = g * p;
    auto phi = euler(p, g);
    auto e = keyGeneration(try_count, halfp.to_long(), false);//getRandom(0, n);
    while (Euclid(e, n) != 1) {
        i++;
        e = keyGeneration(try_count, i, false);//getRandom(0, n);
    }
    auto y = extendedEuclid(phi, e);
    auto d = y % phi;
    auto code = encryption(word, e, n);
    auto res = decryption(code, d, n);
    auto ttime = difftime(time(nullptr), start);
    std::cout << "result for " << ttime << " sec and word: \' " << word << " \' = " << res
              << " is they same? => " << bool(res == word);
    if(p < LONG_MAX && g < LONG_MAX)
        writeSimples(p.to_long(), g.to_long());
    writeTimes(ttime, n, e , d);
    return 0;
}
