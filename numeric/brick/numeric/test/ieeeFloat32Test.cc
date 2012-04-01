/**
***************************************************************************
* @file IEEEFloat32Test.cc
* 
* Source file defining IEEEFloat32Test class.
*
* Copyright (C) 2004-2005,2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <cmath>
#include <map>

#include <brick/numeric/ieeeFloat32.hh>
#include <brick/test/testFixture.hh>

namespace brick {

  namespace numeric {

    class IEEEFloat32Test : public brick::test::TestFixture<IEEEFloat32Test> {

    public:

      IEEEFloat32Test();
      ~IEEEFloat32Test() {};

      void setUp(const std::string& /* testName */) {}
      void tearDown(const std::string& /* testName */) {}

      // Tests of member functions.
      void testConstructor__void();
      void testConstructor__float();
      void testConstructor__uchar_uchar_uchar_uchar();
      void testConstructor__IEEEFloat32();
      void testDestructor();
      void testConversionOperator_NativeFloatType();
      void testGetByte__size_t();
      void testSetValue__float();
      void testSetValue__uchar_uchar_uchar_uchar();

    private:

      /** 
       * This private function sumply updates the internal state to make
       * sure the specified value is one of the ones explicitly tested.
       * 
       * @param value This argument 
       */
      void
      addTestFloat(float floatValue,
                   unsigned char byte0,
                   unsigned char byte1,
                   unsigned char byte2,
                   unsigned char byte3);

    
      // This vector holds an array of floats for which we know the
      // correct IEEE 32 bit representation.
      std::vector<float> m_knownFloatVector;
    
      // These maps are a convenient way of looking up the appropriate
      // binary representation for the floats in m_knownFloatVector,
      // above.
      std::map<float, unsigned char> m_byte0Map;
      std::map<float, unsigned char> m_byte1Map;
      std::map<float, unsigned char> m_byte2Map;
      std::map<float, unsigned char> m_byte3Map;
    
    }; // class IEEEFloat32Test


    /* ============== Member Function Definititions ============== */

    IEEEFloat32Test::
    IEEEFloat32Test()
      : TestFixture<IEEEFloat32Test>("IEEEFloat32Test"),
        m_knownFloatVector(),
        m_byte0Map(),
        m_byte1Map(),
        m_byte2Map(),
        m_byte3Map()
    {
      // Register all tests.
      BRICK_TEST_REGISTER_MEMBER(testConstructor__void);
      BRICK_TEST_REGISTER_MEMBER(testConstructor__float);
      BRICK_TEST_REGISTER_MEMBER(testConstructor__uchar_uchar_uchar_uchar);
      BRICK_TEST_REGISTER_MEMBER(testConstructor__IEEEFloat32);
      BRICK_TEST_REGISTER_MEMBER(testDestructor);
      BRICK_TEST_REGISTER_MEMBER(testConversionOperator_NativeFloatType);
      BRICK_TEST_REGISTER_MEMBER(testGetByte__size_t);
      BRICK_TEST_REGISTER_MEMBER(testSetValue__float);
      BRICK_TEST_REGISTER_MEMBER(testSetValue__uchar_uchar_uchar_uchar);


      // Fill in the ground truth for later tests.
      this->addTestFloat(0.0, 0x00, 0x00, 0x00, 0x00);
      this->addTestFloat(1.0, 0x3f, 0x80, 0x00, 0x00);
      this->addTestFloat(2.0, 0x40, 0x00, 0x00, 0x00);
      this->addTestFloat(static_cast<float>(23.594), 0x41, 0xbc, 0xc0, 0x83);
      this->addTestFloat(static_cast<float>(-123.941), 0xc2, 0xf7, 0xe1, 0xcb);
      this->addTestFloat(1.75, 0x3f, 0xe0, 0x00, 0x00);
      this->addTestFloat(static_cast<float>(-34.432175), 0xc2, 0x09, 0xba, 0x8c);
      this->addTestFloat(0.0625, 0x3d, 0x80, 0x00, 0x00);
    }


    void
    IEEEFloat32Test::
    testConstructor__void()
    {
      IEEEFloat32 ieeeFloat;
      BRICK_TEST_ASSERT(static_cast<float>(ieeeFloat) == 0.0);
      BRICK_TEST_ASSERT(ieeeFloat.getByte(0) == m_byte0Map[0.0]);
      BRICK_TEST_ASSERT(ieeeFloat.getByte(1) == m_byte1Map[0.0]);
      BRICK_TEST_ASSERT(ieeeFloat.getByte(2) == m_byte2Map[0.0]);
      BRICK_TEST_ASSERT(ieeeFloat.getByte(3) == m_byte3Map[0.0]);
    }

    void
    IEEEFloat32Test::
    testConstructor__float()
    {
      for(size_t index0 = 0; index0 < m_knownFloatVector.size(); ++index0) {
        float floatValue = m_knownFloatVector[index0];
        IEEEFloat32 ieeeFloat(floatValue);
        BRICK_TEST_ASSERT(static_cast<float>(ieeeFloat) == floatValue);
        BRICK_TEST_ASSERT(ieeeFloat.getByte(0) == m_byte0Map[floatValue]);
        BRICK_TEST_ASSERT(ieeeFloat.getByte(1) == m_byte1Map[floatValue]);
        BRICK_TEST_ASSERT(ieeeFloat.getByte(2) == m_byte2Map[floatValue]);
        BRICK_TEST_ASSERT(ieeeFloat.getByte(3) == m_byte3Map[floatValue]);
      }
    }
  
    void
    IEEEFloat32Test::
    testConstructor__uchar_uchar_uchar_uchar()
    {
      for(size_t index0 = 0; index0 < m_knownFloatVector.size(); ++index0) {
        float floatValue = m_knownFloatVector[index0];
        unsigned char byte0 = m_byte0Map[floatValue];
        unsigned char byte1 = m_byte1Map[floatValue];
        unsigned char byte2 = m_byte2Map[floatValue];
        unsigned char byte3 = m_byte3Map[floatValue];
        IEEEFloat32 ieeeFloat(byte0, byte1, byte2, byte3);
        BRICK_TEST_ASSERT(static_cast<float>(ieeeFloat) == floatValue);
        BRICK_TEST_ASSERT(ieeeFloat.getByte(0) == byte0);
        BRICK_TEST_ASSERT(ieeeFloat.getByte(1) == byte1);
        BRICK_TEST_ASSERT(ieeeFloat.getByte(2) == byte2);
        BRICK_TEST_ASSERT(ieeeFloat.getByte(3) == byte3);
      }
    }
  
    void
    IEEEFloat32Test::
    testConstructor__IEEEFloat32()
    {
      for(size_t index0 = 0; index0 < m_knownFloatVector.size(); ++index0) {
        float floatValue = m_knownFloatVector[index0];
        IEEEFloat32 ieeeFloat0(floatValue);
        IEEEFloat32 ieeeFloat1(ieeeFloat0);
        BRICK_TEST_ASSERT(static_cast<float>(ieeeFloat1) == floatValue);
        BRICK_TEST_ASSERT(ieeeFloat1.getByte(0) == m_byte0Map[floatValue]);
        BRICK_TEST_ASSERT(ieeeFloat1.getByte(1) == m_byte1Map[floatValue]);
        BRICK_TEST_ASSERT(ieeeFloat1.getByte(2) == m_byte2Map[floatValue]);
        BRICK_TEST_ASSERT(ieeeFloat1.getByte(3) == m_byte3Map[floatValue]);
      }
    }
  
    void
    IEEEFloat32Test::
    testDestructor()
    {
      // Empty.
    }
  
    void
    IEEEFloat32Test::
    testConversionOperator_NativeFloatType()
    {
      // No explicit test.
    }
  
    void
    IEEEFloat32Test::
    testGetByte__size_t()
    {
      // No explicit test.
    }
  
    void
    IEEEFloat32Test::
    testSetValue__float()
    {
      for(size_t index0 = 0; index0 < m_knownFloatVector.size(); ++index0) {
        float floatValue = m_knownFloatVector[index0];
        IEEEFloat32 ieeeFloat;
        ieeeFloat.setValue(floatValue);
        BRICK_TEST_ASSERT(static_cast<float>(ieeeFloat) == floatValue);
        BRICK_TEST_ASSERT(ieeeFloat.getByte(0) == m_byte0Map[floatValue]);
        BRICK_TEST_ASSERT(ieeeFloat.getByte(1) == m_byte1Map[floatValue]);
        BRICK_TEST_ASSERT(ieeeFloat.getByte(2) == m_byte2Map[floatValue]);
        BRICK_TEST_ASSERT(ieeeFloat.getByte(3) == m_byte3Map[floatValue]);
      }
    }
  
    void
    IEEEFloat32Test::
    testSetValue__uchar_uchar_uchar_uchar()
    {
      for(size_t index0 = 0; index0 < m_knownFloatVector.size(); ++index0) {
        float floatValue = m_knownFloatVector[index0];
        unsigned char byte0 = m_byte0Map[floatValue];
        unsigned char byte1 = m_byte1Map[floatValue];
        unsigned char byte2 = m_byte2Map[floatValue];
        unsigned char byte3 = m_byte3Map[floatValue];
        IEEEFloat32 ieeeFloat;
        ieeeFloat.setValue(byte0, byte1, byte2, byte3);
        BRICK_TEST_ASSERT(static_cast<float>(ieeeFloat) == floatValue);
        BRICK_TEST_ASSERT(ieeeFloat.getByte(0) == byte0);
        BRICK_TEST_ASSERT(ieeeFloat.getByte(1) == byte1);
        BRICK_TEST_ASSERT(ieeeFloat.getByte(2) == byte2);
        BRICK_TEST_ASSERT(ieeeFloat.getByte(3) == byte3);
      }
    }


    // ================= Private member functions below =================

    // This private function sumply updates the internal state to make
    // sure the specified value is one of the ones explicitly tested.
    void
    IEEEFloat32Test::
    addTestFloat(float floatValue,
                 unsigned char byte0,
                 unsigned char byte1,
                 unsigned char byte2,
                 unsigned char byte3)
    {
      m_knownFloatVector.push_back(floatValue);
      m_byte0Map[floatValue] = byte0;
      m_byte1Map[floatValue] = byte1;
      m_byte2Map[floatValue] = byte2;
      m_byte3Map[floatValue] = byte3;
    }
  
  } // namespace numeric

} // namespace brick


#if 0

int main(int argc, char** argv)
{
  brick::numeric::IEEEFloat32Test currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::numeric::IEEEFloat32Test currentTest;

}

#endif
