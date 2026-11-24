namespace UM {
    struct PolyLineAttributes {
        PolyLineAttributes() = default;
        PolyLineAttributes(PolyLineAttributes& p)       = default;
        PolyLineAttributes(PolyLineAttributes&& p)      = default;
        PolyLineAttributes(const PolyLineAttributes& p) = default;
        PolyLineAttributes& operator=(const PolyLineAttributes& p) = default;

        PolyLineAttributes(AttributeMap points, AttributeMap edges) : points(std::move(points)), edges(std::move(edges)) {}

        PolyLineAttributes(std::initializer_list<NamedAttribute> list) {
            for (auto& na : list) {
                switch (na.attribute.kind()) {
                    case AttributeBase::POINTS: points[na.name] = na.attribute.get_ptr(); break;
                    case AttributeBase::EDGES:   edges[na.name] = na.attribute.get_ptr(); break;
                    default:
                        um_assert(false);
                }
            }
        }
        AttributeMap points = {};
        AttributeMap edges  = {};
    };

    template <typename T> struct PolyLine::Attribute : GenericAttribute<T> {
        Attribute(T def = T());
        Attribute(PolyLine &seg, T def = T());
        Attribute(const PolyLine &seg, T def = T());
        Attribute(std::string name, PolyLineAttributes &attributes, PolyLine &seg, T def = T());
        bool bind(std::string name, PolyLineAttributes &attributes, PolyLine &seg); // bind on the attribute if found in the collection (return value true), otherwise pushes a new attribute in the collection
        virtual AttributeBase::TYPE kind() const { return AttributeBase::EDGES; }
    };



    template <typename T> using EdgeAttribute = PolyLine::Attribute<T>;
    template <typename T> PointSet::Attribute<T>::Attribute(PolyLine &m, T def) : Attribute(m.points, def) {}
    template <typename T> PointSet::Attribute<T>::Attribute(const PolyLine &m, T def) : Attribute(m.points, def) {}

    template <typename T> PolyLine::Attribute<T>::Attribute(T def) : GenericAttribute<T>(def) {}
    template <typename T> PolyLine::Attribute<T>::Attribute(PolyLine &m, T def) : GenericAttribute<T>(def, m.nedges()) {
        m.attr.push_back(this->ptr);
    }

    template <typename T> PolyLine::Attribute<T>::Attribute(const PolyLine &m, T def) : GenericAttribute<T>(def, m.nedges()) {
    }

    template <typename T> PolyLine::Attribute<T>::Attribute(std::string name, PolyLineAttributes &attributes, PolyLine &m, T def) : GenericAttribute<T>(def) {
        bind_attribute(this, name, m.nedges(), attributes.edges, m.attr);
    }

    template <typename T> bool PolyLine::Attribute<T>::bind(std::string name, PolyLineAttributes &attributes, PolyLine &m) {
        um_assert(!this->bound());
        return bind_attribute(this, name, m.nedges(), attributes.edges, m.attr);
    }

    template <typename T> PointSet::Attribute<T>::Attribute(std::string name, PolyLineAttributes &attributes, PolyLine &seg, T def) : GenericAttribute<T>(def) {
        bind_attribute(this, name, seg.nverts(), attributes.points, *seg.points.attr);
    }

    template <typename T> bool PointSet::Attribute<T>::bind(std::string name, PolyLineAttributes &attributes, PolyLine &seg) {
        um_assert(!this->bound());
        return bind_attribute(this, name, seg.nverts(), attributes.points, *seg.points.attr);
    }
}

