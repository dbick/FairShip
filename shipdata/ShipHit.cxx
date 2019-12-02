#include "ShipHit.h"


// -----   Default constructor   -------------------------------------------
ShipHit::ShipHit()
  : TObject(),
    fdigi(0),
    fDetectorID(-1)
{
}
// -------------------------------------------------------------------------



// -----   Standard constructor   ------------------------------------------
ShipHit::ShipHit(Int_t detID, Float_t digi)
  :TObject(),
   fdigi        (digi),
   fDetectorID  (detID)
{
}

//copy constructor by Stefan
ShipHit::ShipHit(const ShipHit& original)
: TObject(original), fdigi(original.fdigi), fDetectorID(original.fDetectorID)
{

}

// -------------------------------------------------------------------------


// -----   Destructor   ----------------------------------------------------
ShipHit::~ShipHit() { }
// -------------------------------------------------------------------------

//assignment operator by Stefan
ShipHit& ShipHit::operator=(const ShipHit& rhs)
{
	fdigi = rhs.fdigi;
	fDetectorID = rhs.fDetectorID;
	TObject::operator=(rhs);

	return *this;
}

ClassImp(ShipHit)
