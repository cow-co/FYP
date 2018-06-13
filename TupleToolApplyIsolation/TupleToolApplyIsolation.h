#ifndef TupleToolApplyIsolation_H
#define TupleToolApplyIsolation_H 1

// Include files
// from Gaudi
#include "DecayTreeTupleBase/TupleToolBase.h"
#include "Kernel/IParticleTupleTool.h"            // Interface
#include "TrackInterfaces/ITrackVertexer.h"
#include "Kernel/Particle2MCLinker.h"
#include "Kernel/IParticle2MCAssociator.h"
//#include "LoKi/DistanceCalculator.h"
//#include "LoKi/DistanceCalculatorBase.h"
//#include "ChargedParticleMakerBase.h"



class IDVAlgorithm;
class IDistanceCalculator;
class IVertexFit;
class IPVReFitter;

namespace LHCb {
  class Particle;
  class Vertex;
};


/** @class TupleToolVtxIsoln TupleToolVtxIsoln.h
 *
 * \brief Fill isolation information for DecayTreeTuple
 *
 * - head_NOPARTWITHINDCHI2WDW : no. of non-signal particles that when added to vertex give delta chi2 < specified window
 * - head_NOPARTWITHINCHI2WDW : no. of non-signal particles that when added to vertex give chi2 < specified window
 * head_SMALLESTCHI2: chi2 of smallest chi2 combination with any of the input Particles
 * head_SMALLESTDELTACHI2: delta chi2 of smallest delta chi2 combination with any of the input Particles
 *
 * \sa DecayTreeTuple
 *
 *  @todo Maybe one should get Tracks instead of Particles?
 *
 *  @author Mitesh Patel, Patrick Koppenburg
 *  @date   2008-04-15
 */

double MinOfVector(const std::vector<double> vec);

class TupleToolApplyIsolation : public TupleToolBase, virtual public IParticleTupleTool {
public:
  /// Standard constructor
  TupleToolApplyIsolation( const std::string& type,
		     const std::string& name,
		     const IInterface* parent);

  virtual ~TupleToolApplyIsolation( ){}; ///< Destructor

  virtual StatusCode initialize();

  StatusCode fill( const LHCb::Particle*
		   , const LHCb::Particle*
		   , const std::string&
		   , Tuples::Tuple& );

protected:

private:



  Float_t opening,minipchi2,newfdchi2,oldfdchi2,ghostprob,trackchi2,deltafd,pt,ip,chi2,type,vertexchi2,Dst_PT,dummy;

  bool isTrackInDecay(const LHCb::Track*, std::vector<const LHCb::Track*>);
  double getminipchi(const LHCb::Particle*);
  double getfdchi2(const LHCb::Track*, LHCb::Vertex);
  double getopening(const LHCb::Track*,const  LHCb::Particle*);
  const LHCb::Vertex* originVertex( const  LHCb::Particle*, const LHCb::Particle* ) const;

  /**
   * Checks if the particle is a pion with pT < 300 MeV
   * @param  particle The candidate in question
   * @return          1 if positive slow pion, 0 if not a slow pion, -1 if negative slow pion.
   */
  int SlowPionCheck(const LHCb::Particle* particle) const;

  double VertexChi2(const LHCb::Particle& partA, const LHCb::Particle& partB) const;

  /**
   * Calculates the Delta Vertex Chi2 of the candidate with each of the given particles.
   */
  double DVChi2(const LHCb::Particle& candidate, const LHCb::Particle& partA, const LHCb::Particle& partB) const;
  double GetDistanceBetweenPoints(const Gaudi::XYZPoint& pointA, const Gaudi::XYZPoint& pointB, Gaudi::SymMatrix3x3 covarPosA, 
                                  Gaudi::SymMatrix3x3 covarPosB, double& outUncertainty) const;
  double AngleBetweenTracks(const LHCb::Particle& partA, const LHCb::Particle& partB) const;

  IDVAlgorithm* m_dva;

  IDistanceCalculator *m_dist;
  const IVertexFit* m_pVertexFit;
  IParticle2MCAssociator *m_p2mcAssoc;
  IPVReFitter* m_pvReFitter;
  double m_deltaChi2;
  double m_Chi2;
  std::string m_typeVertexFit;
  std::string m_outputSuffix;
  std::string m_weightsName;
  std::vector<std::string> m_inputParticles;
};
#endif // TupleToolApplyIsolation_H
