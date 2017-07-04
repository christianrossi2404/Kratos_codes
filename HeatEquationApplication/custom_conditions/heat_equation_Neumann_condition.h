//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Christian Rossi
//

#ifndef KRATOS_HEAT_EQUATION_NEUMANN_CONDITION
#define KRATOS_HEAT_EQUATION_NEUMANN_CONDITION

// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/condition.h"
#include "includes/model_part.h"
#include "includes/serializer.h"
#include "includes/process_info.h"

// Application includes

#include "includes/deprecated_variables.h"
#include "includes/cfd_variables.h"

namespace Kratos
{
///@addtogroup HeatEquationApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{


template< unsigned int TDim, unsigned int TNumNodes = TDim >
class KRATOS_API(HEAT_EQUATION_APPLICATION) HeatEquationNeumannCondition : public Condition
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of HeatEquationNeumannCondition
    KRATOS_CLASS_POINTER_DEFINITION(HeatEquationNeumannCondition);

    struct ConditionDataStruct
    {
        double wGauss;                  // Gauss point weight
        // double charVel;                 // Problem characteristic velocity (used in the outlet inflow prevention)
        double delta;                   // Non-dimensional positive sufficiently small constant (used in the outlet inflow prevention)
        array_1d<double, 3> Normal;     // Condition normal
        array_1d<double, TNumNodes> N;  // Gauss point shape functions values
    };

    typedef Node < 3 > NodeType;

    typedef Properties PropertiesType;

    typedef Geometry<NodeType> GeometryType;

    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;

    typedef Vector VectorType;

    typedef Matrix MatrixType;

    typedef std::size_t IndexType;

    typedef std::vector<std::size_t> EquationIdVectorType;

    typedef std::vector< Dof<double>::Pointer > DofsVectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    /** Admits an Id as a parameter.
      @param NewId Index for the new         // Struct to pass around the data
        ElementDataStruct data;
        this->FillElementData(data, rCurrentProcessInfo);condition
      */
    HeatEquationNeumannCondition(IndexType NewId = 0):Condition(NewId)
    {
    }

    /// Constructor using an array of nodes
    /**
     @param NewId Index of the new condition
     @param ThisNodes An array containing the nodes of the new condition
     */
    HeatEquationNeumannCondition(IndexType NewId, const NodesArrayType& ThisNodes):
        Condition(NewId,ThisNodes)
    {
    }

    /// Constructor using Geometry
    /**
     @param NewId Index of the new condition
     @param pGeometry Pointer to a geometry object
     */
    HeatEquationNeumannCondition(IndexType NewId, GeometryType::Pointer pGeometry):
        Condition(NewId,pGeometry)
    {
    }

    /// Constructor using Properties
    /**
     @param NewId Index of the new element
     @param pGeometry Pointer to a geometry object
     @param pProperties Pointer to the element's properties
     */
    HeatEquationNeumannCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties):
        Condition(NewId,pGeometry,pProperties)
    {
    }

    /// Copy constructor.
    HeatEquationNeumannCondition(HeatEquationNeumannCondition const& rOther):
        Condition(rOther)
    {
    }

    /// Destructor.
    virtual ~HeatEquationNeumannCondition() {}


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator
    HeatEquationNeumannCondition & operator=(HeatEquationNeumannCondition const& rOther)
    {
        Condition::operator=(rOther);
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    /// Create a new HeatEquationNeumannCondition object.
    /**
      @param NewId Index of the new condition
      @param ThisNodes An array containing the nodes of the new condition
      @param pProperties Pointer to the element's properties
      */
    virtual Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
    {
        return Condition::Pointer(new HeatEquationNeumannCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
    }

    /// Create a new HeatEquationNeumannCondition object.
    /**
      @param NewId Index of the new condition
      @param pGeom A pointer to the condition's geometry
      @param pProperties Pointer to the element's properties
      */
    virtual Condition::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
    {
        return boost::make_shared< HeatEquationNeumannCondition >(NewId, pGeom, pProperties);
    }

    /**
     * Clones the selected element variables, creating a new one
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @return a Pointer to the new element
     */
    virtual Condition::Pointer Clone(IndexType NewId, NodesArrayType const& rThisNodes) const
    {
        Condition::Pointer pNewCondition = Create(NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

        pNewCondition->SetData(this->GetData());
        pNewCondition->SetFlags(this->GetFlags());

        return pNewCondition;
    }

    /// Calculates the LHS and RHS condition contributions
    /**
     * Clones the selected element variables, creating a new one
     * @param rLeftHandSideMatrix: reference to the LHS matrix
     * @param rRightHandSideVector: reference to the RHS matrix
     * @param rCurrentProcessInfo: reference to the ProcessInfo (unused)
     */
    virtual void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                      VectorType& rRightHandSideVector,
                                      ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        constexpr unsigned int MatrixSize = TNumNodes;

        if (rLeftHandSideMatrix.size1() != MatrixSize)
            rLeftHandSideMatrix.resize(MatrixSize, MatrixSize, false); //false says not to preserve existing storage!!

        if (rRightHandSideVector.size() != MatrixSize)
            rRightHandSideVector.resize(MatrixSize, false); //false says not to preserve existing storage!!

        // Struct to pass around the data
        ConditionDataStruct data;

        // Allocate memory needed
        array_1d<double, MatrixSize> rhs_gauss;
        bounded_matrix<double, MatrixSize, MatrixSize> lhs_gauss;

        // LHS and RHS contributions initialization
        noalias(rLeftHandSideMatrix) = ZeroMatrix(MatrixSize,MatrixSize);
        noalias(rRightHandSideVector) = ZeroVector(MatrixSize);

        // Compute condition unit normal vector
        this->CalculateNormal(data.Normal); //this already contains the area
        const double A = norm_2(data.Normal);
        data.Normal /= A;

        // Store the outlet inflow prevention constants in the data structure
        data.delta = 1e-2; // TODO: Decide if this constant should be fixed or not
        // data.charVel = rCurrentProcessInfo[CHARACTERISTIC_VELOCITY];

        // Gauss point information
        GeometryType& rGeom = this->GetGeometry();
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(GeometryData::GI_GAUSS_2);
        const unsigned int NumGauss = IntegrationPoints.size();
        Vector GaussPtsJDet = ZeroVector(NumGauss);
        rGeom.DeterminantOfJacobian(GaussPtsJDet, GeometryData::GI_GAUSS_2);
        const MatrixType Ncontainer = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);

        // Loop on gauss points
        for(unsigned int igauss = 0; igauss<NumGauss; igauss++)
        {
            data.N = row(Ncontainer, igauss);
            const double J = GaussPtsJDet[igauss];
            data.wGauss = J * IntegrationPoints[igauss].Weight();

            ComputeGaussPointRHSContribution(rhs_gauss, data);
            ComputeGaussPointLHSContribution(lhs_gauss, data);

            noalias(rLeftHandSideMatrix) += lhs_gauss;
            noalias(rRightHandSideVector) += rhs_gauss;
        }

        KRATOS_CATCH("")
    }


    virtual void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                       ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        constexpr unsigned int MatrixSize = TNumNodes;

        if (rLeftHandSideMatrix.size1() != MatrixSize)
            rLeftHandSideMatrix.resize(MatrixSize, MatrixSize, false); //false says not to preserve existing storage!!

        // LHS contributions initialization
        noalias(rLeftHandSideMatrix) = ZeroMatrix(MatrixSize,MatrixSize);

        KRATOS_CATCH("")
    }

    /// Calculates the RHS condition contributions
    /**
     * Clones the selected element variables, creating a new one
     * @param rRightHandSideVector: reference to the RHS matrix
     * @param rCurrentProcessInfo: reference to the ProcessInfo (unused)
     */
    virtual void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                        ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        constexpr unsigned int MatrixSize = TNumNodes;

        if (rRightHandSideVector.size() != MatrixSize)
            rRightHandSideVector.resize(MatrixSize, false); //false says not to preserve existing storage!!

        // Struct to pass around the data
        ConditionDataStruct data;

        // Allocate memory needed
        array_1d<double,MatrixSize> rhs_gauss;

        // Loop on gauss points
        noalias(rRightHandSideVector) = ZeroVector(MatrixSize);

        // Compute condition normal
        this->CalculateNormal(data.Normal); //this already contains the area
        const double A = norm_2(data.Normal);
        data.Normal /= A;

        // Store the outlet inflow prevention constants in the data structure
        data.delta = 1e-2; // TODO: Decide if this constant should be fixed or not
        // data.charVel = rCurrentProcessInfo[CHARACTERISTIC_VELOCITY];

        // Gauss point information
        GeometryType& rGeom = this->GetGeometry();
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(GeometryData::GI_GAUSS_2);
        const unsigned int NumGauss = IntegrationPoints.size();
        Vector GaussPtsJDet = ZeroVector(NumGauss);
        rGeom.DeterminantOfJacobian(GaussPtsJDet, GeometryData::GI_GAUSS_2);
        const MatrixType Ncontainer = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);

        for(unsigned int igauss = 0; igauss<NumGauss; igauss++)
        {
            data.N = row(Ncontainer, igauss);
            const double J = GaussPtsJDet[igauss];
            data.wGauss = J * IntegrationPoints[igauss].Weight();

            ComputeGaussPointRHSContribution(rhs_gauss, data);

            noalias(rRightHandSideVector) += rhs_gauss;
        }

        KRATOS_CATCH("")
    }


    /// Condition check
    /**
     * @param rCurrentProcessInfo: reference to the ProcessInfo
     */
    virtual int Check(const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        int Check = Condition::Check(rCurrentProcessInfo); // Checks id > 0 and area > 0

        if (Check != 0)
        {
            return Check;
        }
        else
        {
            // Check that all required variables have been registered

            if(TEMPERATURE.Key() == 0)
                KRATOS_ERROR << "TEMPERATURE Key is 0. Check if the application was correctly registered.";


            // Checks on nodes
            // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
            for(unsigned int i=0; i<this->GetGeometry().size(); ++i)
            {

                if(this->GetGeometry()[i].SolutionStepsDataHas(TEMPERATURE) == false)
                    KRATOS_ERROR << "missing TEMPERATURE variable on solution step data for node " << this->GetGeometry()[i].Id();

            }

            return Check;
        }

        KRATOS_CATCH("");
    }

    /// Provides the global indices for each one of this element's local rows.
    /** This determines the elemental equation ID vector for all elemental DOFs
     * @param rResult A vector containing the global Id of each row
     * @param rCurrentProcessInfo the current process info object (unused)
     */
    virtual void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);

    /// Returns a list of the element's Dofs
    /**
     * @param ElementalDofList the list of DOFs
     * @param rCurrentProcessInfo the current process info instance
     */
    virtual void GetDofList(DofsVectorType& ConditionDofList, ProcessInfo& CurrentProcessInfo);

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "HeatEquationNeumannCondition" << TDim << "D";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "HeatEquationNeumannCondition";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

    void CalculateNormal(array_1d<double,3>& An);

    void ComputeGaussPointLHSContribution(bounded_matrix<double,TNumNodes,TNumNodes>& lhs, const ConditionDataStruct& data);
    void ComputeGaussPointRHSContribution(array_1d<double,TNumNodes>& rhs, const ConditionDataStruct& data);

    void ComputeRHSNeumannContribution(array_1d<double,TNumNodes>& rhs, const ConditionDataStruct& data);
    

    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition );
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition );
    }

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{


    ///@}

}; // Class HeatEquationNeumannCondition


///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template< unsigned int TDim, unsigned int TNumNodes >
inline std::istream& operator >> (std::istream& rIStream, HeatEquationNeumannCondition<TDim,TNumNodes>& rThis)
{
    return rIStream;
}

/// output stream function
template< unsigned int TDim, unsigned int TNumNodes >
inline std::ostream& operator << (std::ostream& rOStream, const HeatEquationNeumannCondition<TDim,TNumNodes>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

///@} addtogroup block


}  // namespace Kratos.

#endif // KRATOS_HEAT_EQUATION_NEUMANN_CONDITION
