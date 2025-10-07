#include "WaveFront.hpp"
#include "utils.hpp"
#include <stdexcept>
#include <cmath>

// Constructor
WaveFront::WaveFront(double size, double pixel_size, ray normal, double wavelength, FieldType source, double w0, int l, int p)
    : size(size), pixel_size(pixel_size), normal(normal), wavelength(wavelength), source(source), w0(w0), l(l), p(p)
{
    N = (int)(size / pixel_size);
    Ex.resize(N, std::vector<std::complex<double>>(N, {0.0, 0.0}));
    Ey.resize(N, std::vector<std::complex<double>>(N, {0.0, 0.0}));
    get_LocalFrame();
}

// Getters
double WaveFront::getSize() { return size; }
double WaveFront::getPixelSize() { return pixel_size; }
double WaveFront::getWavelength() { return wavelength; }
ray WaveFront::getNormal() { return normal; }

// Local frame computation
void WaveFront::get_LocalFrame()
{
    vec3 w = normal.dir();
    vec3 temp = (std::abs(w.x()) < 0.9) ? vec3(1.0, 0.0, 0.0) : vec3(0.0, 1.0, 0.0);

    u = unit_vector(cross(w, temp));
    v = unit_vector(cross(w, u));

    if (w.x() > 0.9)
    {
        temp = u;
        u = v;
        v = temp;
    }
}

// Fresnel propagation
void WaveFront::propagate(double z)
{
    if (z == 0.0)
        return;

    normal.propagate(z);

    const double dx = pixel_size;
    const double k_mag = 2 * PI / wavelength;

    fftw_complex *inp = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N * N);
    fftw_complex *out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N * N);

    if (!inp || !out)
    {
        if (inp)
            fftw_free(inp);
        if (out)
            fftw_free(out);
        throw std::bad_alloc();
    }

    auto process_component = [&](const std::vector<std::vector<std::complex<double>>> &A,
                                 std::vector<std::vector<std::complex<double>>> &Aout)
    {
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
            {
                int kidx = idx(i, j);
                inp[kidx][0] = A[i][j].real();
                inp[kidx][1] = A[i][j].imag();
            }

        fftw_plan forward = fftw_plan_dft_2d(N, N, inp, out, FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_execute(forward);
        fftw_destroy_plan(forward);

        for (int u = 0; u < N; ++u)
        {
            int uu = u - N / 2;
            double fx = double(uu) / (N * dx);
            for (int v = 0; v < N; ++v)
            {
                int vv = v - N / 2;
                double fy = double(vv) / (N * dx);
                double H_phase = -PI * wavelength * z * (fx * fx + fy * fy);
                std::complex<double> H = std::polar(1.0, H_phase);

                int kidx = idx(u, v);
                std::complex<double> S(out[kidx][0], out[kidx][1]);
                S *= H;

                out[kidx][0] = S.real();
                out[kidx][1] = S.imag();
            }
        }

        fftw_plan inverse = fftw_plan_dft_2d(N, N, out, inp, FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_execute(inverse);
        fftw_destroy_plan(inverse);

        double norm = 1.0 / double(N * N);
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
            {
                int kidx = idx(i, j);
                std::complex<double> val(inp[kidx][0] * norm, inp[kidx][1] * norm);
                Aout[i][j] = val;
            }
    };

    auto apply_shift = [&](std::vector<std::vector<std::complex<double>>> &A)
    {
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                A[i][j] *= std::complex<double>((((i + j) & 1) ? -1.0 : 1.0), 0.0);
    };

    std::vector<std::vector<std::complex<double>>> Ex_out(N, std::vector<std::complex<double>>(N));
    std::vector<std::vector<std::complex<double>>> Ey_out(N, std::vector<std::complex<double>>(N));

    auto Ex_copy = Ex;
    apply_shift(Ex_copy);
    process_component(Ex_copy, Ex_out);
    apply_shift(Ex_out);

    auto Ey_copy = Ey;
    apply_shift(Ey_copy);
    process_component(Ey_copy, Ey_out);
    apply_shift(Ey_out);

    Ex.swap(Ex_out);
    Ey.swap(Ey_out);

    fftw_free(inp);
    fftw_free(out);
}

// Phase shift
void WaveFront::phaseShift(double phi)
{
    std::complex<double> ph = std::polar(1.0, phi);
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
        {
            Ex[i][j] *= ph;
            Ey[i][j] *= ph;
        }
}

// Scaling
void WaveFront::scale(double factor)
{
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
        {
            Ex[i][j] *= factor;
            Ey[i][j] *= factor;
        }
}

// Intensity
std::vector<std::vector<double>> WaveFront::Intensity() const
{
    std::vector<std::vector<double>> intensity(N, std::vector<double>(N, 0.0));
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            intensity[i][j] = sq(std::norm(Ex[i][j])) + sq(std::norm(Ey[i][j]));
    return intensity;
}

std::vector<std::vector<double>> WaveFront::Phase() const
{
    std::vector<std::vector<double>> phase(N, std::vector<double>(N, 0.0));
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            phase[i][j] = std::arg(Ex[i][j]);
    return phase;
}

void WaveFront::initialize() {
    double k = 2 * PI / wavelength ;
    double half = size / 2.0 ;

    for (int i = 0 ; i < N ; i++) {
        double x = (i - N / 2) * pixel_size ;
        
    }
}

// Operators
WaveFront WaveFront::operator+(const WaveFront &other)
{
    WaveFront C(this->getSize(), this->getPixelSize(), this->getNormal(), this->getWavelength(), this->source, this->w0, this->l, this->p);
    for (int i = 0; i < C.N; i++)
        for (int j = 0; j < C.N; j++)
        {
            C.Ex[i][j] = this->Ex[i][j] + other.Ex[i][j];
            C.Ey[i][j] = this->Ey[i][j] + other.Ey[i][j];
        }
    return C;
}

WaveFront WaveFront::operator-(const WaveFront &other)
{
    WaveFront C(this->getSize(), this->getPixelSize(), this->getNormal(), this->getWavelength(), this->source, this->w0, this->l, this->p);
    for (int i = 0; i < C.N; i++)
        for (int j = 0; j < C.N; j++)
        {
            C.Ex[i][j] = this->Ex[i][j] - other.Ex[i][j];
            C.Ey[i][j] = this->Ey[i][j] - other.Ey[i][j];
        }
    return C;
}

WaveFront &WaveFront::operator+=(const WaveFront &other)
{
    for (int i = 0; i < this->N; i++)
        for (int j = 0; j < this->N; j++)
        {
            this->Ex[i][j] += other.Ex[i][j];
            this->Ey[i][j] += other.Ey[i][j];
        }
    return *this;
}

WaveFront &WaveFront::operator-=(const WaveFront &other)
{
    for (int i = 0; i < this->N; i++)
        for (int j = 0; j < this->N; j++)
        {
            this->Ex[i][j] -= other.Ex[i][j];
            this->Ey[i][j] -= other.Ey[i][j];
        }
    return *this;
}

// Reflection
void WaveFront::reflect(vec3 n)
{
    normal.reflect(n);
}
