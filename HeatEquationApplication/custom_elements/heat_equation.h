//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Christian Rossi
//

#if !defined(KRATOS_HEAT_EQUATION)
#define  KRATOS_HEAT_EQUATION

// System includes

// External includes
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/serializer.h"
#include "utilities/geometry_utilities.h"
// #include "includes/cfd_variables.h"
// #include "fluid_dynamics_application_variables.h"

// Application includes

namespace Kratos
{

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

// TODO: UPDATE THIS INFORMATION
/**
*/
template< unsigned int TDim, unsigned int TNumNodes = TDim + 1 >
class HeatEquation : public Element
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of
    KRATOS_CLASS_POINTER_DEFINITION(HeatEquation);

    struct ElementDataStruct
    {
        bounded_matrix<double, TNumNodes, TDim> v;
        array_1d<double,TNumNodes> temp, tempn, tempnn, Q;

        bounded_matrix<double, TNumNodes, TDim > DN_DX;
        array_1d<double, TNumNodes > N;
	      double rho;
        double cp;
        double k;
        double bdf0;
        double bdf1;
        double bdf2;
        double h;             // Element size
        double volume;        // In 2D: element area. In 3D: element volume
        double delta_t;       // Only, needed if the temporal dependent term is considered in the subscales
        double dyn_tau_coeff; // Only, needed if the temporal dependent term is considered in the subscales
    };

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.

    HeatEquation(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
    {}

    HeatEquation(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
    {}

    /// Destructor.
    virtual ~HeatEquation() {};


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    Element::Pointer Create(IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY
        return boost::make_shared< HeatEquation < TDim, TNumNodes > >(NewId, this->GetGeometry().Create(rThisNodes), pProperties);
        KRATOS_CATCH("");
    }

    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY
        return boost::make_shared< HeatEquation < TDim, TNumNodes > >(NewId, pGeom, pProperties);
        KRATOS_CATCH("");
    }


    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                              VectorType& rRightHandSideVector,
                              ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        std::cout << "Inside CalculateLocalSystem()" << std::endl;

        constexpr unsigned int MatrixSize = TNumNodes;

        if (rLeftHandSideMatrix.size1() != MatrixSize)
            rLeftHandSideMatrix.resize(MatrixSize, MatrixSize, false); //false says not to preserve existing storage!!

        if (rRightHandSideVector.size() != MatrixSize)
            rRightHandSideVector.resize(MatrixSize, false); //false says not to preserve existing storage!!


        // Struct to pass around the data
        ElementDataStruct data;
        this->FillElementData(data, rCurrentProcessInfo);

        // Allocate memory needed
        bounded_matrix<double,MatrixSize, MatrixSize> lhs_local;
        array_1d<double,MatrixSize> rhs_local;

        // Loop on gauss points
        noalias(rLeftHandSideMatrix) = ZeroMatrix(MatrixSize,MatrixSize);
        noalias(rRightHandSideVector) = ZeroVector(MatrixSize);

        //KRATOS_WATCH(rLeftHandSideMatrix)
        //KRATOS_WATCH(rRightHandSideVector)

        // Gauss point position
        bounded_matrix<double,TNumNodes, TNumNodes> Ncontainer;
        GetShapeFunctionsOnGauss(Ncontainer);

        //KRATOS_WATCH(Ncontainer)

        for(unsigned int igauss = 0; igauss<Ncontainer.size1(); igauss++)
        {
            noalias(data.N) = row(Ncontainer, igauss);

            //KRATOS_WATCH(data.N)

            ComputeGaussPointLHSContribution(lhs_local, data);
            ComputeGaussPointRHSContribution(rhs_local, data);
            // here we assume that all the weights of the gauss points are the same so we multiply at the end by Volume/n_nodes
            noalias(rLeftHandSideMatrix) += lhs_local;
            noalias(rRightHandSideVector) += rhs_local;

            // KRATOS_WATCH(lhs_local)
            // KRATOS_WATCH(rhs_local)
        }

        KRATOS_WATCH(data.volume)

        rLeftHandSideMatrix  *= data.volume/static_cast<double>(TNumNodes);
        rRightHandSideVector *= data.volume/static_cast<double>(TNumNodes);

        KRATOS_CATCH("")
    }




    void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY
        // std::cout << "Inside CalculateRightHandSide()" << std::endl;

        constexpr unsigned int MatrixSize = TNumNodes;

        if (rRightHandSideVector.size() != MatrixSize)
            rRightHandSideVector.resize(MatrixSize, false); //false says not to preserve existing storage!!

        // Struct to pass around the data
        ElementDataStruct data;
        this->FillElementData(data, rCurrentProcessInfo);

        // Allocate memory needed
        array_1d<double,MatrixSize> rhs_local;

        // Gauss point position
        bounded_matrix<double,TNumNodes, TNumNodes> Ncontainer;
        GetShapeFunctionsOnGauss(Ncontainer);

        // Loop on gauss point
        noalias(rRightHandSideVector) = ZeroVector(MatrixSize);
        for(unsigned int igauss = 0; igauss<Ncontainer.size2(); igauss++)
        {
            noalias(data.N) = row(Ncontainer, igauss);

            ComputeGaussPointRHSContribution(rhs_local, data);

            // KRATOS_WATCH(rhs_local)
            //here we assume that all the weights of the gauss points are the same so we multiply at the end by Volume/n_nodes
            noalias(rRightHandSideVector) += rhs_local;
        }


        rRightHandSideVector *= data.volume/static_cast<double>(TNumNodes);

        KRATOS_CATCH("")

    }


    /// Checks the input and that all required Kratos variables have been registered.
    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo The ProcessInfo of the ModelPart that contains this element.
     * @return 0 if no errors were found.
     */
    virtual int Check(const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        // Perform basic element checks
        int ErrorCode = Kratos::Element::Check(rCurrentProcessInfo);
        if(ErrorCode != 0) return ErrorCode;

        // Check that all required variables have been registered
        //if(VELOCITY.Key() == 0)
        //    KRATOS_THROW_ERROR(std::invalid_argument,"VELOCITY Key is 0. Check if the application was correctly registered.","");
        if(TAU.Key() == 0)
           KRATOS_THROW_ERROR(std::invalid_argument,"TAU Key is 0. Check if the application was correctly registered.","");
        if(DENSITY.Key() == 0)
           KRATOS_THROW_ERROR(std::invalid_argument,"DENSITY Key is 0. Check if the application was correctly registered.","");
        if(CONDUCTIVITY.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"CONDUCTIVITY Key is 0. Check if the application was correctly registered.","");
        if(SPECIFIC_HEAT.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"SPECIFIC_HEAT Key is 0. Check if the application was correctly registered.","");

        // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
        for(unsigned int i=0; i<this->GetGeometry().size(); ++i)
        {
            //if(this->GetGeometry()[i].SolutionStepsDataHas(VELOCITY) == false)
            //    KRATOS_THROW_ERROR(std::invalid_argument,"Missing VELOCITY variable on solution step data for node ",this->GetGeometry()[i].Id());
            if(this->GetGeometry()[i].SolutionStepsDataHas(TEMPERATURE) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,"Missing TEMPERATURE variable on solution step data for node ",this->GetGeometry()[i].Id());

            if(this->GetGeometry()[i].HasDofFor(TEMPERATURE) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,"Missing TEMPERATURE component degree of freedom on node ",this->GetGeometry()[i].Id());
        }

        // Check constitutive law
      //  if(mpConstitutiveLaw == nullptr)
      //      KRATOS_ERROR << "The constitutive law was not set. Cannot proceed. Call the navier_stokes.h Initialize() method needs to be called.";

      //  mpConstitutiveLaw->Check(GetProperties(), this->GetGeometry(), rCurrentProcessInfo);

        return 0;

        KRATOS_CATCH("");
    }


      virtual std::string Info() const override
    {
        return "HeatEquation #";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info() << Id();
    }

    /// Print object's data.
    // virtual void PrintData(std::ostream& rOStream) const override

    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static member variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{

    // Constitutive law pointer
    // ConstitutiveLaw::Pointer mpConstitutiveLaw = nullptr;

    // Symbolic function implementing the element   GetDofList<TDim,TNumNodes>
    void GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& rCurrentProcessInfo) override;
    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) override;

    void ComputeGaussPointLHSContribution(bounded_matrix<double,TNumNodes,TNumNodes>& lhs, const ElementDataStruct& data);
    void ComputeGaussPointRHSContribution(array_1d<double,TNumNodes>& rhs, const ElementDataStruct& data);

    //double SubscaleErrorEstimate(const ElementDataStruct& data);

    ///@}
    ///@name Protected Operators
    ///@{

    HeatEquation() : Element()
    {
    }


    // Auxiliar function to fill the element data structure
    void FillElementData(ElementDataStruct& rData, const ProcessInfo& rCurrentProcessInfo)
    {
        // Getting data for the given geometry
        // double Volume; // In 2D cases Volume variable contains the element area

        // std::cout << "Inside FillElementData()..." << std::endl;
        GeometryUtils::CalculateGeometryData(this->GetGeometry(), rData.DN_DX, rData.N, rData.volume);
        // std::cout << "After geometry data" << std::endl;
        // Compute element size
        rData.h = ComputeH(rData.DN_DX);

        //KRATOS_WATCH(rData.h)
        // Database access to all of the variables needed
        const Vector& BDFVector = rCurrentProcessInfo[BDF_COEFFICIENTS];

        rData.delta_t = rCurrentProcessInfo[DELTA_TIME];         // Only, needed if the temporal dependent term is considered in the subscales
        // KRATOS_WATCH(rData.delta_t)
        // KRATOS_WATCH(BDFVector) // TODO: Look at this
        rData.bdf0 = BDFVector[0];
        rData.bdf1 = BDFVector[1];
        rData.bdf2 = BDFVector[2];

        //KRATOS_WATCH(rData.bdf0)
        // rData.dyn_tau_coeff = this->GetProperties()[TAU];  // Only, needed if the temporal dependent term is considered in the subscales
        rData.dyn_tau_coeff = rCurrentProcessInfo[TAU];
        // KRATOS_WATCH(rData.dyn_tau_coeff)
        rData.cp = this->GetProperties()[SPECIFIC_HEAT];
        rData.rho = this->GetProperties()[DENSITY];
        rData.k = this->GetProperties()[CONDUCTIVITY];

        // std::cout << "Before loop" << std::endl;
        for (unsigned int i = 0; i < TNumNodes; i++)
        {

            const array_1d<double,3>& vel = this->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);

           for(unsigned int k=0; k<TDim; k++)
           {
               rData.v(i,k)   = vel[k];

           }

           rData.Q[i] = this->GetGeometry()[i].FastGetSolutionStepValue(HEAT_FLUX);
           rData.temp[i] = this->GetGeometry()[i].FastGetSolutionStepValue(TEMPERATURE);
           rData.tempn[i] = this->GetGeometry()[i].FastGetSolutionStepValue(TEMPERATURE,1);
           rData.tempnn[i] = this->GetGeometry()[i].FastGetSolutionStepValue(TEMPERATURE,2);
        }
        // std::cout << "OK filldata" << std::endl;
    }

    //~ template< unsigned int TDim, unsigned int TNumNodes=TDim+1>
    double ComputeH(boost::numeric::ublas::bounded_matrix<double,TNumNodes, TDim>& DN_DX)
    {
        double h=0.0;
        for(unsigned int i=0; i<TNumNodes; i++)
        {
            double h_inv = 0.0;
            for(unsigned int k=0; k<TDim; k++)
            {
                h_inv += DN_DX(i,k)*DN_DX(i,k);
            }
            h += 1.0/h_inv;
        }
        h = sqrt(h)/static_cast<double>(TNumNodes);
        return h;
    }

    // 3D tetrahedra shape functions values at Gauss points
    void GetShapeFunctionsOnGauss(boost::numeric::ublas::bounded_matrix<double,4,4>& Ncontainer)
    {
        Ncontainer(0,0) = 0.58541020; Ncontainer(0,1) = 0.13819660; Ncontainer(0,2) = 0.13819660; Ncontainer(0,3) = 0.13819660;
        Ncontainer(1,0) = 0.13819660; Ncontainer(1,1) = 0.58541020; Ncontainer(1,2) = 0.13819660; Ncontainer(1,3) = 0.13819660;
        Ncontainer(2,0) = 0.13819660; Ncontainer(2,1) = 0.13819660; Ncontainer(2,2) = 0.58541020; Ncontainer(2,3) = 0.13819660;
        Ncontainer(3,0) = 0.13819660; Ncontainer(3,1) = 0.13819660; Ncontainer(3,2) = 0.13819660; Ncontainer(3,3) = 0.58541020;
    }

    // 2D triangle shape functions values at Gauss points
    void GetShapeFunctionsOnGauss(boost::numeric::ublas::bounded_matrix<double,3,3>& Ncontainer)
    {
        const double one_sixt = 1.0/6.0;
        const double two_third = 2.0/3.0;
        Ncontainer(0,0) = one_sixt; Ncontainer(0,1) = one_sixt; Ncontainer(0,2) = two_third;
        Ncontainer(1,0) = one_sixt; Ncontainer(1,1) = two_third; Ncontainer(1,2) = one_sixt;
        Ncontainer(2,0) = two_third; Ncontainer(2,1) = one_sixt; Ncontainer(2,2) = one_sixt;
    }

    // 3D tetrahedra shape functions values at centered Gauss point
    void GetShapeFunctionsOnUniqueGauss(boost::numeric::ublas::bounded_matrix<double,1,4>& Ncontainer)
    {
        Ncontainer(0,0) = 0.25; Ncontainer(0,1) = 0.25; Ncontainer(0,2) = 0.25; Ncontainer(0,3) = 0.25;
    }

    // 2D triangle shape functions values at centered Gauss point
    void GetShapeFunctionsOnUniqueGauss(boost::numeric::ublas::bounded_matrix<double,1,3>& Ncontainer)
    {
        Ncontainer(0,0) = 1.0/3.0; Ncontainer(0,1) = 1.0/3.0; Ncontainer(0,2) = 1.0/3.0;
    }

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

    virtual void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
    }

    virtual void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
    }

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

};

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
                                    Fluid2DASGS& rThis);
 */
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
                                    const Fluid2DASGS& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }*/
///@}

} // namespace Kratos.

#endif //
