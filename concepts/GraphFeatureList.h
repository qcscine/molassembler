namespace MoleculeManip {

class GraphFeatureList {
private:
  std::vector<
    std::shared_ptr<GraphFeature>
  > _features;

public:
  void add(const std::shared_ptr<GraphFeature>& feature_ptr);

  void remove_involving(const AtomIndexType& a);

  /* Information */
  std::vector<
    std::shared_ptr<GraphFeature>
  > get_all_matching_type(
    const std::string& type
  ) const;

  /* Operators */
  std::vector<
    std::shared_ptr<GraphFeature>
  >::const_iterator begin() const {
    return _features.begin();
  }
  std::vector<
    std::shared_ptr<GraphFeature>
  >::const_iterator end() const {
    return _features.end();
  }
};

} // eo namespace
