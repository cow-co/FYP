// Include files
// from Gaudi
#include "GaudiKernel/ToolFactory.h"

// local
#include "TupleToolApplyIsolation.h"
#include <Kernel/GetIDVAlgorithm.h>
#include <Kernel/IDVAlgorithm.h>
#include <Kernel/IDistanceCalculator.h>
#include <Kernel/IVertexFit.h>
#include "GaudiAlg/Tuple.h"
#include "GaudiAlg/TupleObj.h"
#include "Event/Particle.h"
#include "Event/MCParticle.h"
#include "GaudiKernel/DeclareFactoryEntries.h"
#include <functional>
#include "TrackInterfaces/ITrackVertexer.h"
#include "Linker/LinkerTable.h"
#include "Kernel/IPVReFitter.h"

#include <map>
#include <string>
#include <cfloat>   //For DBL_MAX
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TPluginManager.h"
#include "TStopwatch.h"

// Declaration of the Tool Factory
// actually acts as a using namespace TupleTool
DECLARE_TOOL_FACTORY( TupleToolApplyIsolation );
using namespace LHCb;

double MinOfVector(const std::vector<double> vec)
{
    double min = vec[0];

    for(unsigned int i = 1; i < vec.size(); ++i)
    {
        if(vec[i] < min)
        {
            min = vec[i];
        }
    }

    return min;
}

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
TupleToolApplyIsolation::TupleToolApplyIsolation( const std::string& type,
                                      const std::string& name,
                                      const IInterface* parent )
  : TupleToolBase ( type, name , parent )
    , m_dva(0)
    , m_dist(0)
    , m_pVertexFit(0)
{
  declareInterface<IParticleTupleTool>(this);
  m_inputParticles.push_back("/Event/Phys/StdAllNoPIDsPions");
  m_inputParticles.push_back("/Event/Phys/StdNoPIDsUpPions");
  m_inputParticles.push_back("Phys/StdNoPIDsVeloPions");
  //m_inputParticles.push_back("/Event/Phys/StdNoPIDsVeloElectrons");
  //havent removed / added any of this yet
  declareProperty( "MaxDeltaChi2", m_deltaChi2 = 9.0);
  declareProperty( "MaxChi2", m_Chi2 = 9.0);
  declareProperty( "VertexFit", m_typeVertexFit = "default");
  declareProperty("InputParticles", m_inputParticles );
  declareProperty("OutputSuffix", m_outputSuffix = "" );

}
//=============================================================================
StatusCode TupleToolApplyIsolation::initialize() {
  if( ! TupleToolBase::initialize() ) return StatusCode::FAILURE;
  m_dva = Gaudi::Utils::getIDVAlgorithm ( contextSvc() ) ;
  if (0==m_dva) return Error("Couldn't get parent DVAlgorithm",
                             StatusCode::FAILURE);
  m_dist       = tool<IDistanceCalculator>("LoKi::DistanceCalculator",this);
  m_p2mcAssoc  = tool<IParticle2MCAssociator>("DaVinciSmartAssociator", this);
  if( !m_dist ){
    Error("Unable to retrieve the IDistanceCalculator tool");
    return StatusCode::FAILURE;
  }
  m_pvReFitter  = tool<IPVReFitter>("AdaptivePVReFitter", this );
  m_pVertexFit  = tool<IVertexFit>("LoKi::VertexFitter",this);
  //m_pVertexFit= m_dva->vertexFitter();
  //m_pVertexFit= tool<ITrackVertexer>
  if( !m_pVertexFit ){
    Error("Unable to retrieve the IVertexFit tool");
    return StatusCode::FAILURE;
  }
  return StatusCode::SUCCESS;
}
//=============================================================================
StatusCode TupleToolApplyIsolation::fill( const Particle* mother
                                        , const Particle* P
                                        , const std::string& head
                                        , Tuples::Tuple& tuple )
{
    const std::string prefix=fullName(head);
    Assert(P && mother && m_dist,
        "This should not happen, you are inside TupleToolApplyIsoln.cpp :(");
    bool test = true;

    std::vector<const LHCb::Track*> daughtertracks;
    daughtertracks.clear();
    LHCb::Particle::ConstVector source; //representing basically the B0.
                                        //The original particle.
    LHCb::Particle::ConstVector target; //TODO: What is this representing?
    LHCb::Particle::ConstVector finalStates;    //LHCb::Particle::ConstVector is
                                                //a typedef around
                                                //std::vector<const Particle*>
    LHCb::Particle::ConstVector parts2Vertex;
    LHCb::Particle::ConstVector parts2VertexD;

    vertexchi2 = P->endVertex()->chi2();
    parts2Vertex.clear();
    parts2VertexD.clear();

    //

    /**
     * If the particle has no daughters.
     * This if/else checks that the B0 is in the top of the decay chain.
     */
    if (P->isBasicParticle())
    {
        source.push_back(mother);
    }

    else
    {
        source.push_back(P);
    }

    LHCb::Vertex dv2;

    do
    {
        target.clear();

        /**
         * Iterate over the source particles (the B0s)
         */
        for(LHCb::Particle::ConstVector::const_iterator isource = source.begin();
            isource != source.end(); isource++)
        {
            /**
             * If the particle has daughters...
             */
            if(!((*isource)->daughters().empty()))
            {
                //...Get the daughters, and then...
                LHCb::Particle::ConstVector tmp = (*isource)->daughtersVector();

                /**
                 * ...Iterate over the daughters, and for each one...
                 */
                for(LHCb::Particle::ConstVector::const_iterator itmp = tmp.begin();
                    itmp!=tmp.end(); itmp++)
                {
                    target.push_back(*itmp);    //...Add to the target list, before...

                    /**
                     * ...Adding the final states, i.e. particles with proto and
                     * ignoring photons.  Protoparticle is the precursor class
                     * to Particle (stores some raw detector data).  PDG code
                     * (i.e. particleID().pid()) 22 is the code for a photon:
                     * http://pdg.lbl.gov/2006/reviews/montecarlorpp.pdf
                     */
                    if((*itmp)->proto() && 22 != (*itmp)->particleID().pid())
                    {
                        /**
                         * Add the particle to the final states and its track to
                         * the daughter tracks
                         */
                        finalStates.push_back(*itmp);
                        daughtertracks.push_back((*itmp)->proto()->track());
                    }
                }
            } // if endVertex
        } // isource

        source = target;
    }

    while(target.size() > 0);   //Each time we run the do loop, the target vector
                                //gets smaller, because it is being copied to the
                                //source vector, and the source vector determines
                                //what we get in the target, so it shrinks.

    if (msgLevel(MSG::DEBUG))
    {
        debug() << "Final states size= " <<  finalStates.size()  << endreq;
    }

    // now push back particles of decay chain into a vector
    // to (re)-perform the vertex fit.
    LHCb::Vertex v;

    if(P->isBasicParticle())
    {
        parts2Vertex.push_back(P);
    }

    /**
     * If the particle has daughters, do a vertex fit with those.
     */
    else
    {
        parts2Vertex = P->daughtersVector();
        StatusCode sc = m_pVertexFit->fit(v, parts2Vertex);
    }

    LHCb::Particle::ConstVector theParts;
    LHCb::Particle muPlus, muMinus; //we will always have two objects tagged as muons (at least,
                                    //with the toolset we're using).

    /**
     * Grab the muons so we can vertex them with our other tracks
     */
    for(LHCb::Particle::ConstVector::iterator it = finalStates.begin();
        it != finalStates.end(); ++it)
    {
        int pdgID = (*it)->particleID().pid();

        /**
         * If the PDG ID matches that of a muon, add it to the vector.
         */
        if(pdgID == 13)
        {
            muMinus = (**it);
        }

        else if(pdgID == -13)
        {
            muPlus = (**it);
        }
    }

    //Our discriminating variables
    unsigned int numToSave = 7; //how many candidates should we save?  This determines the number of
                                //branches per variable.
    std::vector<double> minipchi2s(numToSave, 0.); //fills with zeroes

    std::vector<double> maxPTs(numToSave, DBL_MIN);  //The highest track p_T values.  These tracks are the ones for which we'll take the DVChi2
    std::vector<double> deltaVertexChi2s(numToSave, 0.);

    std::vector<double> fdchi2s(numToSave, 0.);
    std::vector<double> minDVChi2Chi2s(2, 0.);  //the values of the individual Chi2s contributing to the minimum DVCHi2.

    double minDVChi2 = DBL_MAX;

    Gaudi::XYZPoint posOfMinChi2[2]; //positions of the vertices of the particle with the lowest Chi2 with the negative muon, with each muon.
    double rMmX1MpX1 = -1.0; //magnitude of the position vector of the vertex with the lowest Chi2 with the - muon
    double uncertaintyRMmX1MpX1 = -1.0;

    LHCb::Vertex muPlusVtx, muMinusVtx; //vertices for each muon.
    double minChi2Plus = DBL_MAX, minChi2Minus = DBL_MAX; //min Chi2s with each muon.
    Gaudi::XYZPoint posOfMinChi2Plus, posOfMinChi2Minus;  //positions of essentially the vertices each muon comes from (i.e. tracks that fit best with each mu)
    double uncertaintyMpPos = -1.0, uncertaintyMmPos = -1.0;  //variance on the muon vertex positions.
    
    double distFromPosVtxToMuVtx = -1.0, uncertaintyDistMpMM = -1.0; //distance from posOfMinChi2Plus to the dimuon vertex, and the variance on it
    double distFromNegVtxToMuVtx = -1.0, uncertaintyDistMmMM = -1.0; //distance from posOfMinChi2Minus to the dimuon vertex, and its variance
    
    double ipMmX1ToMp = -1.0; //the point of closest approach on the + muon's path to the - muon/X1 vertex.

    double distBetweenMinMuVerts = -1.0, uncertaintyDistMX1MX2 = -1.0; //distance between the minimum verts formed by each muon.
    
    LHCb::Vertex tripleVertex;  //the vertex used to construct vertices from three (rather than the usual two) tracks.
    double dimuonX1Chi2 = -1.0; //the vertex chi2 of the vertex formed by the two muons and the track with the best vertex with the - muon.
    double dimuonX1Dist = -1.0; //the distance of the abovementioned vertex from the origin of the coordinate system. 
    double dimuonX1Uncertainty = -1.0;
    double dimuonX1DistScaled = -1.0;  //the above distance, but scaled by the uncertainty on the vertex position. 
    double deltaDimuonPosWithX1 = -1.0; //the change in vertex position when we add in the extra track to the dimuon vertex
    double deltaDimuonPosWithX1Scaled = -1.0; //the change in vertex position when we add in the extra track to the dimuon vertex

    LHCb::Vertex muonVertex;  //dimuon vertex.
    double dimuonDist = -1.0, dimuonDistScaled = -1.0, uncertaintyDimuonDist = -1.0;
    
    Gaudi::SymMatrix3x3 dimuonVtxCovMatrix;
    Gaudi::SymMatrix3x3 muX1VtxCovMatrix;    
    Gaudi::SymMatrix3x3 muX2VtxCovMatrix;    
    Gaudi::SymMatrix3x3 dimuonX1VtxCovMatrix;
    Gaudi::SymMatrix3x3 blankMatrix;
    
    Gaudi::XYZPoint origin(0.0, 0.0, 0.0);  //invent an origin position, with 0 covariances (for use in the distance function we have)
    
    m_pVertexFit->fit(muonVertex, muMinus, muPlus); //fit the vertex
    dimuonVtxCovMatrix = muonVertex.covMatrix();
    dimuonDist = GetDistanceBetweenPoints(origin, muonVertex.position(), blankMatrix, dimuonVtxCovMatrix, uncertaintyDimuonDist);
    dimuonDistScaled = dimuonDist / uncertaintyDimuonDist;
    
    double muonVtxChi2 = ((VertexBase)muonVertex).chi2PerDoF(); //vertex quality
    double ipDimuon = -1.0;
    m_dist->distance(&muPlus, &muMinus, ipDimuon);

    /**
     * Iterate over the input particle strings, and for each one...
     * FIXME: Could optimise this by minimising the number of variables initialised on each iteration.
     */
    for(std::vector<std::string>::iterator i = m_inputParticles.begin();
        i !=m_inputParticles.end(); ++i)
    {
        /**
         * ...Check that there are, in fact, particles in the file(?), and then...
         */
        if(!exist<LHCb::Particle::Range>(*i+"/Particles"))
        {
            if (msgLevel(MSG::DEBUG))
            {
                debug() << "No particles at " << *i << " !1!1!" << endreq;
            }

            continue;
        }

        //...Grab those particles, before...
        LHCb::Particle::Range parts = get<LHCb::Particle::Range>(*i+"/Particles");

        if (msgLevel(MSG::DEBUG))
        {
            debug() << "Getting particles from " << *i
                    << " with " << (parts).size() << " particles" << endreq;
        }

        /**
         * ...Iterating over them, allowing us to...
         */
        for(LHCb::Particle::Range::const_iterator iparts = (parts).begin();
            iparts != (parts).end(); ++iparts)
        {
            const LHCb::Particle* part = (*iparts);

            /**
             * ...only pick particular types of track that are not in the decay
             * chain...
             * (5 denotes a Downstream track, and recall that daughtertracks was
             * filled with the final-state tracks).
             */
            if(part->proto()->track()->type() < 5 &&
                !isTrackInDecay(part->proto()->track(), daughtertracks))
            {
                LHCb::Vertex vtxWithExtraTrack;

                // add track to the vertex formed by
                // particles of the original decay chain and re-perform
                // the vertex fit
                parts2Vertex.push_back(*iparts);
                StatusCode sc3 = m_pVertexFit->fit(vtxWithExtraTrack,parts2Vertex);

                // after the fit remove that track
                parts2Vertex.pop_back();

                // now calculate a bunch of variables related to the extra track
                // and the particles of the decay chain in question
                opening   = getopening(part->proto()->track(), P);
                minipchi2 = getminipchi(part);
                newfdchi2 = getfdchi2(part->proto()->track(), vtxWithExtraTrack);
                oldfdchi2 = getfdchi2(part->proto()->track(), v);
                ghostprob = part->proto()->track()->ghostProbability();
                trackchi2 = part->proto()->track()->chi2PerDoF();
                deltafd   = log10(fabs(newfdchi2 - oldfdchi2)) - 7;
                type      = part->proto()->track()->type();

                if (newfdchi2 - oldfdchi2 < 0)
                {
                    deltafd = deltafd * -1.;
                }

                newfdchi2 = log10(newfdchi2);

                /**
                 * Now we apply some constraints
                 */

                // VeLo only tracks dont have momentum information so treat
                // accordingly
                if(part->proto()->track()->type() == 1)
                {
                    pt = part->proto()->track()->momentum().z();
                }

                else
                {
                    pt = part->proto()->track()->pt();
                }

                // only look at tracks not coming from primary vertex
                if(minipchi2 < 5)
                {
                    continue;
                }

                // only look at tracks of good quality
                if(trackchi2 > 3)
                {
                    continue;
                }

                // make sure track is not made up of wrong hits and also
                // not a clone of the original track
                if(ghostprob > 0.5)
                {
                    continue;
                }

                if(part->proto()->track()->type() == 3 && !(opening > 0.994 ))
                {
                    continue;
                }

                if(part->proto()->track()->type() == 4 && !(opening > 0.98))
                {
                    continue;
                }

                if(part->proto()->track()->type() == 1 && !(opening > 0.98))
                {
                    continue;
                }

                StatusCode sc = StatusCode::SUCCESS;
                double tmpip, tmpchi2;
                StatusCode dump = m_dist->distance((const LHCb::Particle *) part,
                                                    (const LHCb::Vertex *)&v,
                                                    tmpip,tmpchi2);
                chi2 = tmpchi2;

                double chi2Minus = VertexChi2(*part, muMinus);
                double chi2Plus = VertexChi2(*part, muPlus);
		
                if(chi2Plus == 0.0 || chi2Minus == 0.0)
                {
                  chi2Plus = -1.0;
                  chi2Minus = -1.0;
                }
                
                if(chi2Plus > -1.0 && chi2Minus > -1.0)
                {
                  double dvChi2 = DVChi2(*part, muMinus, muPlus);

                  //Put the tracks into various lists, sorted in descending order of track pt.
                  for(unsigned int i = 0; i < numToSave; ++i)
                  {
                    if(pt > maxPTs[i])
                    {
                      //Scooch all the lower values down by one place.
                      for(unsigned int j = numToSave - 1; j > i; --j)
                      {
                        maxPTs[j] = maxPTs[j - 1];
                        deltaVertexChi2s[j] = deltaVertexChi2s[j - 1];
                        minipchi2s[j] = minipchi2s[j - 1];
                      }
                      
                      maxPTs[i] = pt;
                      deltaVertexChi2s[i] = dvChi2;
                      minipchi2s[i] = minipchi2;

                      break;
                    }
                  }

                  if(dvChi2 < minDVChi2 && dvChi2 > 0.0)
                  {
                    minDVChi2 = dvChi2;
                    minDVChi2Chi2s[0] = chi2Minus;
                    minDVChi2Chi2s[1] = chi2Plus;
                  }

                  //Find the particles with the lowest chi2 with each muon.
                  if(chi2Minus < minChi2Minus && chi2Minus > 0.0)
                  {
                    LHCb::Vertex posVtx;  //the vertex of this track with the positive muon.  In BG, this is NOT the same as the 
                                          //from which the positive muon originates.

                    minChi2Minus = chi2Minus;

                    //Vertex our track with each muon and find the distance between those vertices.
                    m_pVertexFit->fit(muMinusVtx, *part, muMinus);
                    m_pVertexFit->fit(posVtx, *part, muPlus);

                    posOfMinChi2[0] = muMinusVtx.position();
                    posOfMinChi2[1] = posVtx.position();

                    rMmX1MpX1 = GetDistanceBetweenPoints(posOfMinChi2[0], posOfMinChi2[1], muMinusVtx.covMatrix(), posVtx.covMatrix(), uncertaintyRMmX1MpX1);

                    posOfMinChi2Minus = muMinusVtx.position();
                    muX1VtxCovMatrix = muMinusVtx.covMatrix();
                    
                    m_dist->distance(part, (const VertexBase*)&muMinusVtx, ipMmX1ToMp);
                    
                    //semi-arbitrarily, we choose this track as the one which we add in to the muon vertex.
                    StatusCode sc = m_pVertexFit->fit(tripleVertex, *part, muPlus, muMinus);
                    
                    if(sc != StatusCode::SUCCESS)
                    {
                      if(msgLevel(MSG::DEBUG))
                      {
                        debug() << "Triple-Track Vertex Fit Failed!" << endreq;
                      }
                     
                      dimuonX1Chi2 = -1.0;
                      dimuonX1Dist = -1.0;
                      dimuonX1DistScaled = -1.0;
                    }
                    
                    else
                    {
                      dimuonX1Chi2 = ((VertexBase)tripleVertex).chi2PerDoF();
                      dimuonX1Dist = GetDistanceBetweenPoints(origin, tripleVertex.position(), blankMatrix,
                                                              tripleVertex.covMatrix(), dimuonX1Uncertainty);

                      double delWithXUncertainty = -1.0;                      
                      deltaDimuonPosWithX1 = GetDistanceBetweenPoints(tripleVertex.position(), muonVertex.position(), tripleVertex.covMatrix(), 
                                                                      muonVertex.covMatrix(), delWithXUncertainty);
                      dimuonX1VtxCovMatrix = tripleVertex.covMatrix();
                      dimuonX1DistScaled = dimuonX1Dist / dimuonX1Uncertainty;
                      deltaDimuonPosWithX1Scaled = deltaDimuonPosWithX1 / delWithXUncertainty;
                    }
                  }

                  if(chi2Plus < minChi2Plus && chi2Plus > 0.0)
                  {
                    minChi2Plus = chi2Plus;
                    m_pVertexFit->fit(muPlusVtx, *part, muPlus);
                    posOfMinChi2Plus = muPlusVtx.position();
                    muX2VtxCovMatrix = muPlusVtx.covMatrix();
                  }
                }
            }
        } // end particles loop
    }//end particle types loop
    
    //Only calculate and output the variables if there were no errors in the loop.
    if(minChi2Minus > -1.0 && minChi2Plus > -1.0 && minDVChi2 > -1.0)
    {
      if(dimuonX1Chi2 > -1.0 && dimuonX1Dist > -1.0)
      {
        tuple->column(prefix + "_MMX1_CHI2", dimuonX1Chi2);
        tuple->column(prefix + "_DELTA_CHI2_MMX1_MM", fabs(dimuonX1Chi2 - muonVtxChi2));
        tuple->column(prefix + "_R_MMX1", dimuonX1Dist);
        tuple->column(prefix + "_R_MMX1_SCALED", dimuonX1DistScaled);
        tuple->column(prefix + "_DEL_R_MM_MMX1", deltaDimuonPosWithX1);
        tuple->column(prefix + "_DEL_R_MM_MMX1_SCALED", deltaDimuonPosWithX1Scaled);
      }
      
      if(ipMmX1ToMp != -1.0)
      {
        tuple->column(prefix + "_IP_MmX1_Mp", ipMmX1ToMp);
      }
      
      if(ipDimuon != -1.0)
      {
        tuple->column(prefix + "_IP_DIMUON", ipDimuon);
      }
      
      distFromPosVtxToMuVtx = GetDistanceBetweenPoints(posOfMinChi2Plus, muonVertex.position(), muPlusVtx.covMatrix(), dimuonVtxCovMatrix, uncertaintyDistMpMM);
      distFromNegVtxToMuVtx = GetDistanceBetweenPoints(posOfMinChi2Minus, muonVertex.position(), muMinusVtx.covMatrix(), dimuonVtxCovMatrix, uncertaintyDistMmMM);

      distBetweenMinMuVerts = GetDistanceBetweenPoints(posOfMinChi2Minus, posOfMinChi2Plus, muMinusVtx.covMatrix(), muPlusVtx.covMatrix(), uncertaintyDistMX1MX2);

      if(distBetweenMinMuVerts >= 0.0 &&
         distFromPosVtxToMuVtx >= 0.0 &&
         distFromNegVtxToMuVtx >= 0.0 &&
         rMmX1MpX1             >= 0.0)
      {
        tuple->column(prefix + "_R_MpX1_MmX2", distBetweenMinMuVerts);
        tuple->column(prefix + "_R_MM_MpX", distFromPosVtxToMuVtx);
        tuple->column(prefix + "_R_MpX_MmX", rMmX1MpX1);
        
        tuple->column(prefix + "_R_MpX1_MmX2_X", fabs(posOfMinChi2Minus.x() - posOfMinChi2Plus.x()));
        tuple->column(prefix + "_R_MpX1_MmX2_Y", fabs(posOfMinChi2Minus.y() - posOfMinChi2Plus.y()));
        tuple->column(prefix + "_R_MpX1_MmX2_Z", fabs(posOfMinChi2Minus.z() - posOfMinChi2Plus.z()));
        
        tuple->column(prefix + "_R_MM_MpX_X", fabs(posOfMinChi2Plus.x() - muonVertex.position().x()));
        tuple->column(prefix + "_R_MM_MpX_Y", fabs(posOfMinChi2Plus.y() - muonVertex.position().y()));
        tuple->column(prefix + "_R_MM_MpX_Z", fabs(posOfMinChi2Plus.z() - muonVertex.position().z()));
        
        tuple->column(prefix + "_R_MM_MmX_X", fabs(posOfMinChi2Minus.x() - muonVertex.position().x()));
        tuple->column(prefix + "_R_MM_MmX_Y", fabs(posOfMinChi2Minus.y() - muonVertex.position().y()));
        tuple->column(prefix + "_R_MM_MmX_Z", fabs(posOfMinChi2Minus.z() - muonVertex.position().z()));
      }

      for(size_t is = 0; is < numToSave; ++is)
      {
          char branch_count[512];
          sprintf(branch_count,"_%d",is);
          tuple->column(prefix + "_ISOLATION_MINIPCHI2_" + std::string(branch_count) + m_outputSuffix,  minipchi2s[is]);

          if(deltaVertexChi2s[is] >= 0.0)
          {
            tuple->column(prefix + "_DELTA_VERTEX_CHI2_" + std::string(branch_count) + m_outputSuffix,  deltaVertexChi2s[is]);
          }
      }

      if(minDVChi2 >= 0.0)
      {
        tuple->column(prefix + "_MIN_DV_CHI2_" + m_outputSuffix, minDVChi2);
        tuple->column(prefix + "_MIN_DV_CHI2_CHI2_0" + m_outputSuffix, minDVChi2Chi2s[0]);
        tuple->column(prefix + "_MIN_DV_CHI2_CHI2_1" + m_outputSuffix, minDVChi2Chi2s[1]);
      }
      
      distBetweenMinMuVerts /= uncertaintyDistMX1MX2;
      distFromPosVtxToMuVtx /= uncertaintyDistMpMM;
      rMmX1MpX1 /= uncertaintyRMmX1MpX1;

      if(distBetweenMinMuVerts >= 0.0 &&
         distFromPosVtxToMuVtx >= 0.0 &&
         distFromNegVtxToMuVtx >= 0.0 &&
         rMmX1MpX1         >= 0.0)
      {
        tuple->column(prefix + "_R_MpX1_MmX2_SCALED", distBetweenMinMuVerts);
        tuple->column(prefix + "_R_MM_MpX_SCALED", distFromPosVtxToMuVtx);
        tuple->column(prefix + "_R_MpX_MmX_SCALED", rMmX1MpX1);
        
        tuple->matrix(prefix + "_VARIANCE_MM_POS", dimuonVtxCovMatrix, 3, 3);
        tuple->matrix(prefix + "_VARIANCE_MX1_POS", muX1VtxCovMatrix, 3, 3);
        tuple->matrix(prefix + "_VARIANCE_MX2_POS", muX2VtxCovMatrix, 3, 3);
        tuple->matrix(prefix + "_VARIANCE_MMX1_POS", dimuonX1VtxCovMatrix, 3, 3);
      }
    }
    
    //These should always be ok since we will always have two muons AT LEAST.
    tuple->column(prefix + "_MUON_VTX", muonVtxChi2);
    tuple->column(prefix + "_R_MM", dimuonDist);
    tuple->column(prefix + "_R_MM_SCALED", dimuonDistScaled);
    
    return StatusCode(test);
}

double TupleToolApplyIsolation::AngleBetweenTracks(const LHCb::Particle& partA, const LHCb::Particle& partB) const
{
  //Get the momentum 3-vectors.
  Gaudi::XYZPoint momVecA(partA.momentum().Px(), partA.momentum().Py(), partA.momentum().Pz());
  Gaudi::XYZPoint momVecB(partB.momentum().Px(), partB.momentum().Py(), partB.momentum().Pz());

  //Normalise the momenta to get their unit vectors.
  momVecA /= (sqrt(momVecA.mag2()));
  momVecB /= (sqrt(momVecB.mag2()));

  //Dot product of the momenta directions.
  double dot = momVecA.X() * momVecB.X() + momVecA.Y() * momVecB.Y() + momVecA.Z() * momVecB.Z();

  //angle between the momenta.
  return acos(dot);
}

double TupleToolApplyIsolation::GetDistanceBetweenPoints(const Gaudi::XYZPoint& pointA, const Gaudi::XYZPoint& pointB, Gaudi::SymMatrix3x3 covarPosA,
                                                          Gaudi::SymMatrix3x3 covarPosB, double& outUncertainty) const
{
  Gaudi::XYZPoint vecBetween(pointB.X() - pointA.X(), pointB.Y() - pointA.Y(), pointB.Z() - pointA.Z());
  double dist = sqrt(vecBetween.mag2());
  
  double dfByDx = ((pointB.X() - pointA.X()) / dist);
  double dfByDy = ((pointB.Y() - pointA.Y()) / dist);
  double dfByDz = ((pointB.Z() - pointA.Z()) / dist);
  
  //TODO: Could potentially optimise if the covar matrices are symmetrical (which I think they are).  Playing it safe for now though.
  //NOTE: This assumes that the covariances between the different positions are zero (i.e. covar x1x2 = 0 for example).
  double x1x1 = dfByDx * dfByDx * covarPosA[0][0];
  double x1y1 = dfByDx * dfByDy * covarPosA[0][1];
  double x1z1 = dfByDx * dfByDz * covarPosA[0][2];
  
  double y1x1 = dfByDy * dfByDx * covarPosA[1][0];
  double y1y1 = dfByDy * dfByDy * covarPosA[1][1];
  double y1z1 = dfByDy * dfByDz * covarPosA[1][2];
  
  double z1x1 = dfByDz * dfByDx * covarPosA[2][0];
  double z1y1 = dfByDz * dfByDy * covarPosA[2][1];
  double z1z1 = dfByDz * dfByDz * covarPosA[2][2];  
  
  double x2x2 = dfByDx * dfByDx * covarPosB[0][0];
  double x2y2 = dfByDx * dfByDy * covarPosB[0][1];
  double x2z2 = dfByDx * dfByDz * covarPosB[0][2];
  
  double y2x2 = dfByDy * dfByDx * covarPosB[1][0];
  double y2y2 = dfByDy * dfByDy * covarPosB[1][1];
  double y2z2 = dfByDy * dfByDz * covarPosB[1][2];
  
  double z2x2 = dfByDz * dfByDx * covarPosB[2][0];
  double z2y2 = dfByDz * dfByDy * covarPosB[2][1];
  double z2z2 = dfByDz * dfByDz * covarPosB[2][2];
  
  outUncertainty = sqrt(x1x1 + x1y1 + x1z1 + y1x1 + y1y1 + y1z1 + z1x1 + z1y1 + z1z1 + x2x2 + x2y2 + x2z2 + y2x2 + y2y2 + y2z2 + z2x2 + z2y2 + z2z2);

  return dist;
}

double TupleToolApplyIsolation::VertexChi2(const LHCb::Particle& partA, const LHCb::Particle& partB) const
{
    double result = -1.0;
    LHCb::Vertex vtx;
    StatusCode sc = m_pVertexFit->fit(vtx, partA, partB);

    if(sc != StatusCode::SUCCESS)
    {
        if(msgLevel(MSG::DEBUG))
        {
            debug() << "Vertex Fit Failed!" << endreq;
        }
    }

    else
    {
        result = ((VertexBase)vtx).chi2PerDoF();
    }

    return result;
}

double TupleToolApplyIsolation::DVChi2(const LHCb::Particle& candidate, const LHCb::Particle& partA, const LHCb::Particle& partB) const
{
    double delta;
    double chi2PartA = VertexChi2(candidate, partA);
    double chi2PartB = VertexChi2(candidate, partB);

    if(chi2PartA == -1.0 || chi2PartB == -1.0)
    {
        delta = -1.0;
    }

    else
    {
        delta = fabs(chi2PartA - chi2PartB);
    }

    return delta;
}

//=========================================================================
//
//=========================================================================
const Vertex* TupleToolApplyIsolation::originVertex(const Particle* top
                                                    , const Particle* P )
  const {
  if( top == P || P->isBasicParticle() ) return 0;
  const SmartRefVector< LHCb::Particle >& dau = top->daughters ();
  if( dau.empty() ){
    return 0;
  }
  SmartRefVector< LHCb::Particle >::const_iterator it;
  for( it = dau.begin(); dau.end()!=it; ++it ){
    if( P == *it ){ // I found the daughter
      return top->endVertex();
    }
  }
  // vertex not yet found, get deeper in the decay:
  for( it = dau.begin(); dau.end()!=it; ++it ){
    if( P != *it && !(*it)->isBasicParticle() ){
      const Vertex* vv = originVertex( *it, P );
      if( vv ) return vv;
    }
  }
  return 0;
}
//=============================================================================
// Check if the track is already in the decay
//============================================================================
bool TupleToolApplyIsolation::isTrackInDecay(const LHCb::Track* track,
                                             std::vector<const LHCb::Track*>
                                             daughters)
{
    bool isInDecay = false;
    //loop over daughters
    for(std::vector<const LHCb::Track*>::iterator it = daughters.begin();
      it != daughters.end(); ++it)
    {
        const LHCb::Track* itrack = (*it);

        if(itrack)
        {
            if(itrack == track)
            {
                if(msgLevel(MSG::DEBUG))
                {
                    debug() << "Track is in decay, skipping it" << endmsg;
                }

                isInDecay = true;
            }
        }
    }

    return isInDecay;
}
//=============================================================================
// MINIPCHI2 for a track
//=============================================================================
double TupleToolApplyIsolation::getminipchi(const LHCb::Particle* track)
{
    double minchi2 = -1 ;
    const RecVertex::Range PV = m_dva->primaryVertices();   //Get the PVs from the DaVinci algorithm.

    if(!PV.empty())
    {
        /**
         * Iterate over all the reconstructed PVs.
         */
        for(RecVertex::Range::const_iterator pv = PV.begin();
            pv != PV.end(); ++pv)
        {
            double ip, chi2;    //Impact parameter and the change in Chi2 of PV fit when the given track is added to the fit.
            StatusCode test2 = m_dist->distance((const LHCb::Particle*)track, *pv, ip, chi2);   //Fill the IP and Chi2

            /**
             * If the Chi2 is less than the minimum, or if this is the first iteration (in which case
             * minchi2 == -1)...
             */
            if((chi2 < minchi2) || (minchi2 < 0.))
            {
                /**
                 * ...Calculate the new minimum chi2.
                 */
                LHCb::RecVertex newPV(**pv);
                StatusCode scfit = m_pvReFitter->remove(track, &newPV); //remove this track from the PV
                LHCb::RecVertex* newPVPtr = (LHCb::RecVertex*)&newPV;
                test2 = m_dist->distance((LHCb::Particle *)track, (LHCb::VertexBase*)newPVPtr, ip, chi2);
                minchi2 = chi2;
            }
        }
    }

    return minchi2;
}
double TupleToolApplyIsolation::getfdchi2(const LHCb::Track* track,
                                          LHCb::Vertex Vtx){
  double minchi2 = -1 ;
  double fdchi2 = -1;
  double fd;
  const RecVertex::Range PV = m_dva->primaryVertices();
  if ( !PV.empty() ){
    for ( RecVertex::Range::const_iterator pv = PV.begin() ;
          pv!=PV.end() ; ++pv){
      double ip, chi2;
      StatusCode test2 = m_dist->distance ( (const LHCb::Track *)track,
                                            *pv, ip, chi2 );
      if ((chi2<minchi2) || (minchi2<0.))
      {
        minchi2 = chi2 ;
        StatusCode test2 = m_dist->distance ( *pv, &Vtx, fd, fdchi2 );
      }
    }
  }
  return fdchi2;
}
//=============================================================================
// Opening angle for a track and particle
//============================================================================
double TupleToolApplyIsolation::getopening(const LHCb::Track* track,
                                           const  LHCb::Particle* P){
  Gaudi::XYZVector A = P->momentum().Vect();
  Gaudi::XYZVector B = track->momentum();
  double cosopening = A.Dot( B ) / std::sqrt( A.Mag2()*B.Mag2() );
  return cosopening;
}
