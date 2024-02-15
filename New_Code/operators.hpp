/*
  Jan 31, 2024  Jiazheng Sun

  Define ladder operators
*/

#ifndef ORI_SDP_GS_OPERATORS_HPP
#define ORI_SDP_GS_OPERATORS_HPP


template <typename T>
class LadderOp {
protected:
  T idx;
  bool creatorF;
public:
  LadderOp() : idx(idx_invalid<T>()), creatorF(false) {}
  LadderOp(T const & i_idx, bool i_creatorF) : idx(i_idx), creatorF(i_creatorF) {}
}






#endif
