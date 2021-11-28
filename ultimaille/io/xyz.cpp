#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <iterator>
#include "ultimaille/io/xyz.h"


namespace UM {
    void write_xyz(const std::string filename, const PointSet &ps) {
        std::fstream out;
        out.open(filename, std::ios_base::out);
        if (out.fail()) {
            throw std::runtime_error("Failed to open " + filename);
        }
//      out << ps.size() << std::endl;
        out << std::setprecision(std::numeric_limits<double>::max_digits10);
        for (const vec3 &p : ps)
            out << p.x << " " << p.y << " " << p.z << std::endl;
        out.close();
    }

    PointSetAttributes read_xyz(const std::string filename, PointSet &m) {
        m = PointSet();
        std::ifstream in;
        in.open(filename, std::ifstream::in);
        if (in.fail()) {
            throw std::runtime_error("Failed to open " + filename);
        }

        std::string line;
//      int npts = 0;
        bool firstline = true;
        double x,y,z;
        while (!in.eof()) {
            std::getline(in, line);
            std::istringstream iss(line);
            int nfields = std::distance(std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>());
            if (nfields == 1 && firstline) {
                // Some xyz files do have number of points specified, some do not. Ignore that number; push_back into the pointset will do the job just fine.
                // npts = std::stoi(line);
                continue;
            }

            if (nfields>=3) {
                iss.clear();
                iss.seekg(0);
                iss >> x >> y >> z;
                m.push_back({x, y, z});
            }

            if (nfields>6 || (nfields>0 && nfields<3 && !firstline))
                std::cerr << "Warning: suspect number of fields per line (" << nfields << ")" << std::endl;
            firstline = false;
        }

        in.close();
        return {};
    }
}

