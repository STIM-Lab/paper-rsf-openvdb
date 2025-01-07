#include <openvdb/openvdb.h>
#include <openvdb/io/File.h>
#include <openvdb/tools/GridTransformer.h>
#include <iostream>
#include <tira/volume.h>
#include "tira/image/colormap.h"
#include "tira/image.h"
#include <openvdb/tools/VolumeToMesh.h>
#include <openvdb/tools/MeshToVolume.h>
#include <openvdb/tools/Fastsweeping.h>
#include <openvdb/tools/GridOperators.h>
#include <openvdb/tools/LevelsetFilter.h>
#include <openvdb/tools/Filter.h>
#include <openvdb/tools/LevelSetSphere.h>
#include <openvdb/tools/Composite.h>
#include <openvdb/tools/ChangeBackground.h>
#include <limits>
#include <openvdb/tools/Morphology.h> 
#include <openvdb/tools/GridTransformer.h>
#include <openvdb/tools/LevelSetUtil.h>
#include <openvdb/tools/ValueTransformer.h>
#include <cmath>
#include <openvdb/tools/Interpolation.h>
#include <openvdb/tools/Morphology.h> 
#include <mutex>
#include <tbb/enumerable_thread_specific.h>
#include <queue>
#include <tbb/tbb.h>
#include <unordered_set>


using namespace openvdb;



//vdb to image conversion for 3D
template<class GridType> void vdb2img3D(GridType& grid, tira::volume<float>& img) {
    using ValueT = typename GridType::ValueType;

    typename GridType::Accessor accessor = grid.getAccessor();

    //openvdb::Coord dim = grid.evalActiveVoxelDim();
    openvdb::Coord ijk;
    int& i = ijk[0], & j = ijk[1], & k = ijk[2];
    for (i = 0; i < img.X(); i++) {
        for (j = 0; j < img.Y(); j++) {
            for (k = 0; k < img.Z(); k++) {
                float pixel = (float)accessor.getValue(ijk);
                img(i, j, k) = pixel;
            }
        }
    }
}

openvdb::FloatGrid::Ptr heaviside_return(openvdb::FloatGrid::Ptr& grid, float& backgroundvalue)
{
    float epsilon = 0.2;
    // Get accessor for the grid
    openvdb::FloatGrid::Ptr resultGrid = openvdb::FloatGrid::create(backgroundvalue);

    resultGrid->setTransform(grid->transform().copy());

    openvdb::FloatGrid::Accessor resultAccessor = resultGrid->getAccessor();
    // Iterate over the grid's active values
    for (openvdb::FloatGrid::ValueOnIter iter = grid->beginValueOn(); iter; ++iter) {

        float value = *iter;
        // Calculate the Heaviside function value
        float heaviside = (0.5 * (1 + (2 / 3.14159265358979323846) * atan(value / epsilon)));
        //resultGrid->tree().setValueOn(iter.getCoord(), heaviside);
        resultAccessor.setValueOn(iter.getCoord(), heaviside);
    }

    return resultGrid;
}



// Heaviside function
template <class GridType>
void heaviside(GridType& grid)
{
    float epsilon = 0.2;
    // Get accessor for the grid
    typename GridType::Accessor accessor = grid.getAccessor();
    // Iterate over the grid's active values
    for (typename GridType::ValueOnIter iter = grid.beginValueOn(); iter; ++iter) {
        // Get the current value
        const typename GridType::ValueType value = iter.getValue();
        // Calculate the Heaviside function value
        const typename GridType::ValueType heaviside = (0.5 * (1 + (2 / 3.14159265358979323846) * atan(value / epsilon)));
        accessor.setValue(iter.getCoord(), heaviside);
    }
}



// derivative of Heaviside function
template <class GridType>
void deri_heaviside(GridType& grid)
{
    float epsilon = 0.2;

    typename GridType::Accessor accessor = grid.getAccessor();

    for (typename GridType::ValueOnIter iter = grid.beginValueOn(); iter; ++iter) {
        const typename GridType::ValueType value = iter.getValue();
        const typename GridType::ValueType deri_heaviside = (1 / 3.14159265358979323846) * (epsilon / ((epsilon * epsilon) + (value * value)));
        accessor.setValue(iter.getCoord(), deri_heaviside);
    }
}



// inverse a vdb file
template <class GridType>
void vdb_inverse(GridType& grid)
{
    typename GridType::Accessor accessor = grid.getAccessor();

    for (typename GridType::ValueOnIter iter = grid.beginValueOn(); iter; ++iter) {
        const typename GridType::ValueType value = iter.getValue();
        const typename GridType::ValueType inverse_vdb = (1 - value);
        accessor.setValue(iter.getCoord(), inverse_vdb);
    }
}


// Function to multiply two VDB grids
openvdb::FloatGrid::Ptr multiplyGrids(openvdb::FloatGrid::Ptr& grid1, openvdb::FloatGrid::Ptr& grid2, float& backgroundvalue) {
   
    openvdb::FloatGrid::Ptr resultGrid = openvdb::FloatGrid::create(backgroundvalue);

    // Get accessors for the grids
    openvdb::FloatGrid::Accessor resultAccessor = resultGrid->getAccessor();
    openvdb::FloatGrid::Accessor grid1Accessor = grid1->getAccessor();
    openvdb::FloatGrid::Accessor grid2Accessor = grid2->getAccessor();

   
    for (openvdb::FloatGrid::ValueOnIter iter = grid1->beginValueOn(); iter; ++iter) {
        const openvdb::Coord& coord = iter.getCoord();

        // Multiply the values from grid1 and grid2
        resultAccessor.setValue(coord, grid1Accessor.getValue(coord) * grid2Accessor.getValue(coord));
    }

    return resultGrid;
}




// Function to divide two VDB grids
openvdb::FloatGrid::Ptr divideGrids(openvdb::FloatGrid::Ptr& grid1, openvdb::FloatGrid::Ptr& grid2) {
    
    openvdb::FloatGrid::Ptr resultGrid = openvdb::FloatGrid::create();
    resultGrid = grid1->deepCopy();

    // Get accessors for the grids
    openvdb::FloatGrid::Accessor resultAccessor = resultGrid->getAccessor();
    openvdb::FloatGrid::Accessor grid1Accessor = grid1->getAccessor();
    openvdb::FloatGrid::Accessor grid2Accessor = grid2->getAccessor();

   
    for (openvdb::FloatGrid::ValueOnIter iter = grid1->beginValueOn(); iter; ++iter) {
        const openvdb::Coord& coord = iter.getCoord();

        // value from the second grid
        float grid2Value = grid2Accessor.getValue(coord);

        // if the divisor is not zero to avoid division by zero
        if (grid2Value != 0.0f) {
            // Divide the values from grid1 by grid2, and set in resultGrid
            resultAccessor.setValue(coord, grid1Accessor.getValue(coord) / grid2Value);
        }
        else {
            // Handle division by zero 
            resultAccessor.setValue(coord, 0.0f); 
        }
    }

    return resultGrid;
}




openvdb::FloatGrid::Ptr calculateGradientX(const openvdb::FloatGrid::Ptr& grid) {
   /* openvdb::FloatGrid::Ptr gradXGrid = openvdb::FloatGrid::create(0.0);
    gradXGrid->setTransform(grid->transform().copy());*/

    openvdb::FloatGrid::Ptr gradXGrid = openvdb::FloatGrid::create();
    gradXGrid = grid->deepCopy();

    openvdb::FloatGrid::ConstAccessor accessor = grid->getConstAccessor();

    for (openvdb::FloatGrid::ValueOnCIter iter = grid->cbeginValueOn(); iter; ++iter) {
        const openvdb::Coord& xyz = iter.getCoord();
        float valueLeft = accessor.isValueOn(xyz.offsetBy(-1, 0, 0)) ? accessor.getValue(xyz.offsetBy(-1, 0, 0)) : iter.getValue();
        float valueRight = accessor.isValueOn(xyz.offsetBy(1, 0, 0)) ? accessor.getValue(xyz.offsetBy(1, 0, 0)) : iter.getValue();
        float gradX = (valueRight - valueLeft) / 2.0f;
        gradXGrid->tree().setValueOn(xyz, gradX);
    }

    return gradXGrid;
}

openvdb::FloatGrid::Ptr calculateGradientY(const openvdb::FloatGrid::Ptr& grid) {
    /*openvdb::FloatGrid::Ptr gradYGrid = openvdb::FloatGrid::create(0.0);
    gradYGrid->setTransform(grid->transform().copy());*/

    openvdb::FloatGrid::Ptr gradYGrid = openvdb::FloatGrid::create();
    gradYGrid = grid->deepCopy();

    openvdb::FloatGrid::ConstAccessor accessor = grid->getConstAccessor();

    for (openvdb::FloatGrid::ValueOnCIter iter = grid->cbeginValueOn(); iter; ++iter) {
        const openvdb::Coord& xyz = iter.getCoord();
        float valueDown = accessor.isValueOn(xyz.offsetBy(0, -1, 0)) ? accessor.getValue(xyz.offsetBy(0, -1, 0)) : iter.getValue();
        float valueUp = accessor.isValueOn(xyz.offsetBy(0, 1, 0)) ? accessor.getValue(xyz.offsetBy(0, 1, 0)) : iter.getValue();
        float gradY = (valueUp - valueDown) / 2.0f;
        gradYGrid->tree().setValueOn(xyz, gradY);
    }

    return gradYGrid;
}

openvdb::FloatGrid::Ptr calculateGradientZ(const openvdb::FloatGrid::Ptr& grid) {
    /*openvdb::FloatGrid::Ptr gradZGrid = openvdb::FloatGrid::create(0.0);
    gradZGrid->setTransform(grid->transform().copy());*/
    openvdb::FloatGrid::Ptr gradZGrid = openvdb::FloatGrid::create();
    gradZGrid = grid->deepCopy();

    openvdb::FloatGrid::ConstAccessor accessor = grid->getConstAccessor();

    for (openvdb::FloatGrid::ValueOnCIter iter = grid->cbeginValueOn(); iter; ++iter) {
        const openvdb::Coord& xyz = iter.getCoord();
        float valueBelow = accessor.isValueOn(xyz.offsetBy(0, 0, -1)) ? accessor.getValue(xyz.offsetBy(0, 0, -1)) : iter.getValue();
        float valueAbove = accessor.isValueOn(xyz.offsetBy(0, 0, 1)) ? accessor.getValue(xyz.offsetBy(0, 0, 1)) : iter.getValue();
        float gradZ = (valueAbove - valueBelow) / 2.0f;
        gradZGrid->tree().setValueOn(xyz, gradZ);
    }

    return gradZGrid;
}



openvdb::FloatGrid::Ptr calculateGradientMagnitude(
    const openvdb::FloatGrid::Ptr& gradXGrid,
    const openvdb::FloatGrid::Ptr& gradYGrid,
    const openvdb::FloatGrid::Ptr& gradZGrid) {

    /*openvdb::FloatGrid::Ptr gradMagGrid = openvdb::FloatGrid::create(0.0);
    gradMagGrid->setTransform(gradXGrid->transform().copy());*/

    openvdb::FloatGrid::Ptr gradMagGrid = openvdb::FloatGrid::create();
    gradMagGrid = gradXGrid->deepCopy();

    //  gradient grids
    openvdb::FloatGrid::ConstAccessor accessorX = gradXGrid->getConstAccessor();
    openvdb::FloatGrid::ConstAccessor accessorY = gradYGrid->getConstAccessor();
    openvdb::FloatGrid::ConstAccessor accessorZ = gradZGrid->getConstAccessor();

    
    for (openvdb::FloatGrid::ValueOnCIter iter = gradXGrid->cbeginValueOn(); iter; ++iter) {
        const openvdb::Coord& xyz = iter.getCoord();
        // for the current voxel
        float gradX = accessorX.getValue(xyz);
        float gradY = accessorY.getValue(xyz);
        float gradZ = accessorZ.getValue(xyz);
        //  magnitude 
        float magnitude = openvdb::math::Sqrt(gradX * gradX + gradY * gradY + gradZ * gradZ);
        // computed magnitude  in the output grid
        gradMagGrid->tree().setValueOn(xyz, magnitude);
    }

    return gradMagGrid;
}




////////////////////////////////////////TBB_PARALLEL_FOR()///////////////////////////////////////////////////////

// Functor for use with tbb::parallel_for() that operates on a grid's leaf nodes
//TBB_HEAVISIDE
struct Heaviside_Processor {
    using GridType = openvdb::FloatGrid;
    using TreeType = GridType::TreeType;
    using LeafNode = TreeType::LeafNodeType;
    // splits the iteration space of a leaf iterator
    using IterRange = openvdb::tree::IteratorRange<TreeType::LeafCIter>;

    float epsilon;

    Heaviside_Processor(float eps) : epsilon(eps) {}

    void operator()(IterRange& range) const {
        for (; range; ++range) {
            auto& leaf = const_cast<LeafNode&>(*range.iterator()); // modifiable leaf node
            for (auto iter = leaf.beginValueOn(); iter; ++iter) {
                const auto value = iter.getValue();
                const auto heaviside = static_cast<float>(
                    0.5 * (1 + (2 / 3.14159265358979323846) * atan(value / epsilon)));
                leaf.setValueOn(iter.offset(), heaviside);
            }
        }
    }
};

void heaviside_tbb(openvdb::FloatGrid::Ptr grid) {
    Heaviside_Processor proc(0.2); // Epsilon value

    // Wrap a leaf iterator in an IteratorRange and run in parallel
    Heaviside_Processor::IterRange range(grid->tree().cbeginLeaf());
    tbb::parallel_for(range, proc);
}

//Derivative Heavisisde_TBB
struct deriv_Heaviside_Processor {
    using GridType = openvdb::FloatGrid;
    using TreeType = GridType::TreeType;
    using LeafNode = TreeType::LeafNodeType;
    // Isplits the iteration space of a leaf iterator
    using IterRange = openvdb::tree::IteratorRange<TreeType::LeafCIter>;

    float epsilon;

    deriv_Heaviside_Processor(float eps) : epsilon(eps) {}

    void operator()(IterRange& range) const {
        for (; range; ++range) {
            auto& leaf = const_cast<LeafNode&>(*range.iterator()); // modifiable leaf node
            for (auto iter = leaf.beginValueOn(); iter; ++iter) {
                const auto value = iter.getValue();
                const auto deri_heaviside = static_cast<float>(
                    (1 / 3.14159265358979323846) * (epsilon / ((epsilon * epsilon) + (value * value))));
                leaf.setValueOn(iter.offset(), deri_heaviside);
            }
        }
    }
};

void deri_heaviside_tbb(openvdb::FloatGrid::Ptr grid) {
    deriv_Heaviside_Processor proc(0.2); // Epsilon value

    // Wrap a leaf iterator in an IteratorRange and run in parallel
    deriv_Heaviside_Processor::IterRange range(grid->tree().cbeginLeaf());
    tbb::parallel_for(range, proc);
}


// INVERSE VDB
template <class GridType>
struct InverseProcessor {
    using ValueType = typename GridType::ValueType;
    using TreeType = typename GridType::TreeType;
    using LeafNode = typename TreeType::LeafNodeType;
    // splits the iteration space of a leaf iterator
    using IterRange = openvdb::tree::IteratorRange<typename TreeType::LeafCIter>;

    void operator()(IterRange& range) const {
        for (; range; ++range) {
            auto& leaf = const_cast<LeafNode&>(*range.iterator()); // modifiable leaf node
            for (auto iter = leaf.beginValueOn(); iter; ++iter) {
                ValueType value = iter.getValue();
                ValueType inverse_value = 1 - value;  //  inverse of the current value
                leaf.setValueOn(iter.offset(), inverse_value);
            }
        }
    }
};

template <class GridType>
void vdb_inverse_tbb(GridType& grid) {
    InverseProcessor<GridType> proc;

    // Wrap a leaf iterator in an IteratorRange and run in parallel
    typename InverseProcessor<GridType>::IterRange range(grid.tree().cbeginLeaf());
    tbb::parallel_for(range, proc);
}


struct GridMultiplier {
    using GridType = openvdb::FloatGrid;
    using TreeType = GridType::TreeType;
    using LeafNode = TreeType::LeafNodeType;
    using IterRange = openvdb::tree::IteratorRange<typename TreeType::LeafCIter>;

    GridType::ConstAccessor acc1;
    GridType::ConstAccessor acc2;

    GridMultiplier(GridType::ConstAccessor accessor1, GridType::ConstAccessor accessor2)
        : acc1(accessor1), acc2(accessor2) {}

    void operator()(IterRange& range) const {
        for (; range; ++range) {
            auto& leaf1 = const_cast<LeafNode&>(*range.iterator());  
            if (const LeafNode* leaf2 = acc2.probeConstLeaf(leaf1.origin())) {  // Match leaf in grid2
                for (auto iter = leaf1.beginValueOn(); iter; ++iter) {
                    const openvdb::Coord coord = iter.getCoord();
                    if (leaf2->isValueOn(coord)) {  // if the value is active in grid2
                        float value1 = iter.getValue();
                        float value2 = leaf2->getValue(coord);
                        leaf1.setValueOn(iter.offset(), value1 * value2);  // Multiply values
                    }
                    else {
                        leaf1.setValueOff(iter.offset());  // Optionally turn off if no corresponding value in grid2
                    }
                }
            }
            else {
                // no corresponding leaf in grid2
                for (auto iter = leaf1.beginValueOn(); iter; ++iter) {
                    leaf1.setValueOff(iter.offset());
                }
            }
        }
    }
};

openvdb::FloatGrid::Ptr multiplyGrids_tbb(openvdb::FloatGrid::Ptr& grid1, openvdb::FloatGrid::Ptr& grid2, float backgroundValue) {
    openvdb::FloatGrid::Ptr resultGrid = openvdb::FloatGrid::create(backgroundValue);
    
    resultGrid = grid1->deepCopy();
    openvdb::FloatGrid::Accessor resultAccessor = resultGrid->getAccessor();
    GridMultiplier proc(resultGrid->getConstAccessor(), grid2->getConstAccessor());

    GridMultiplier::IterRange range(resultGrid->tree().cbeginLeaf());
    tbb::parallel_for(range, proc);

    return resultGrid;
}


struct VoxelProcessor {
    using GridType = openvdb::FloatGrid;
    using TreeType = GridType::TreeType;
    using LeafNode = TreeType::LeafNodeType;
    using IterRange = openvdb::tree::IteratorRange<typename TreeType::LeafCIter>;

    GridType::ConstPtr input_grid, fout, fin;
    tbb::enumerable_thread_specific<GridType::Ptr>& threadLocalEout;
    tbb::enumerable_thread_specific<GridType::Ptr>& threadLocalEin;
    float factor;

    VoxelProcessor(GridType::ConstPtr ig, GridType::ConstPtr fout, GridType::ConstPtr fin, tbb::enumerable_thread_specific<GridType::Ptr>& eo, tbb::enumerable_thread_specific<GridType::Ptr>& ei, float f)
        : input_grid(ig), fout(fout), fin(fin), threadLocalEout(eo), threadLocalEin(ei), factor(f) {}

    void operator()(IterRange& range) const {
        for (; range; ++range) {
            const LeafNode& leaf = *range.iterator();
            openvdb::Coord xyz;

            auto& EoutGrid = threadLocalEout.local();
            auto& EinGrid = threadLocalEin.local();
            auto EoutAccessor = EoutGrid->getAccessor();
            auto EinAccessor = EinGrid->getAccessor();

            for (auto voxelIter = leaf.beginValueOn(); voxelIter; ++voxelIter) {
                xyz = voxelIter.getCoord();
                float s1 = 0.0f, s2 = 0.0f;
                int sigma1 = 0;  

                for (int u = -sigma1; u <= sigma1; ++u) {
                    for (int v = -sigma1; v <= sigma1; ++v) {
                        for (int w = -sigma1; w <= sigma1; ++w) {
                            openvdb::Coord neighbor = xyz.offsetBy(u, v, w);
                            if (input_grid->tree().isValueOn(neighbor)) {
                                float inputValue = input_grid->tree().getValue(neighbor);
                                float foutValue = 0.0f, finValue = 0.0f;
                                if (fout->tree().isValueOn(neighbor)) {
                                    foutValue = fout->tree().getValue(neighbor);
                                }
                                if (fin->tree().isValueOn(neighbor)) {
                                    finValue = fin->tree().getValue(neighbor);
                                }

                                s1 += factor * std::pow(inputValue - foutValue, 2);
                                s2 += factor * std::pow(inputValue - finValue, 2);
                            }
                        }
                    }
                }

                EoutAccessor.setValueOn(xyz, s1);
                EinAccessor.setValueOn(xyz, s2);
            }
        }
    }
};

void parallelVoxelComputation(openvdb::FloatGrid::ConstPtr input_grid, openvdb::FloatGrid::ConstPtr fout, openvdb::FloatGrid::ConstPtr fin, openvdb::FloatGrid::Ptr Eout, openvdb::FloatGrid::Ptr Ein, float sigma1) {
    float kernelSize = (2 * sigma1) + 1;
    float factor = 1 / std::pow(kernelSize, 3);

    // Create thread-local grids
    tbb::enumerable_thread_specific<openvdb::FloatGrid::Ptr> threadLocalEout([&]() {
        auto localGrid = openvdb::FloatGrid::create();
        localGrid->setTransform(input_grid->transform().copy());
        return localGrid;
        });

    tbb::enumerable_thread_specific<openvdb::FloatGrid::Ptr> threadLocalEin([&]() {
        auto localGrid = openvdb::FloatGrid::create();
        localGrid->setTransform(input_grid->transform().copy());
        return localGrid;
        });

    // iterator range over leaf nodes
    VoxelProcessor::IterRange range(input_grid->tree().cbeginLeaf());

    // tbb 
    tbb::parallel_for(range, VoxelProcessor(input_grid, fout, fin, threadLocalEout, threadLocalEin, factor));

    // merging the thread-local grids 
    for (auto& localGrid : threadLocalEout) {
        Eout->tree().merge(localGrid->tree());
    }

    for (auto& localGrid : threadLocalEin) {
        Ein->tree().merge(localGrid->tree());
    }
}

struct GridDivider {
    using GridType = openvdb::FloatGrid;
    using TreeType = GridType::TreeType;
    using LeafNode = TreeType::LeafNodeType;
    using IterRange = openvdb::tree::IteratorRange<typename TreeType::LeafCIter>;

    GridType::ConstAccessor acc1;
    GridType::ConstAccessor acc2;

    GridDivider(GridType::ConstAccessor accessor1, GridType::ConstAccessor accessor2)
        : acc1(accessor1), acc2(accessor2) {}

    void operator()(IterRange& range) const {
        for (; range; ++range) {
            auto& leaf1 = const_cast<LeafNode&>(*range.iterator());  // modifiable leaf node from grid1
            if (const LeafNode* leaf2 = acc2.probeConstLeaf(leaf1.origin())) {  // matching --- leaf in grid2
                for (auto iter = leaf1.beginValueOn(); iter; ++iter) {
                    const openvdb::Coord coord = iter.getCoord();
                    if (leaf2->isValueOn(coord)) {  // if the value is active in grid2
                        float value1 = iter.getValue();
                        float value2 = leaf2->getValue(coord);
                        if (value2 != 0.0f) {  // no division by zero
                            leaf1.setValueOn(iter.offset(), value1 / value2);  // Divide values
                        }
                        else {
                            leaf1.setValueOff(iter.offset());  // Optionally turn off if divisor is zero
                        }
                    }
                    else {
                        leaf1.setValueOff(iter.offset());  // turn off if no corresponding value in grid2
                    }
                }
            }
            else {
                // if no corresponding leaf in grid2
                for (auto iter = leaf1.beginValueOn(); iter; ++iter) {
                    leaf1.setValueOff(iter.offset());
                }
            }
        }
    }
};

openvdb::FloatGrid::Ptr divideGrids_tbb(openvdb::FloatGrid::Ptr& grid1, openvdb::FloatGrid::Ptr& grid2) {
    openvdb::FloatGrid::Ptr resultGrid = openvdb::FloatGrid::create();
    resultGrid = grid1->deepCopy();
    openvdb::FloatGrid::Accessor resultAccessor = resultGrid->getAccessor();
    GridDivider proc(resultGrid->getConstAccessor(), grid2->getConstAccessor());

    GridDivider::IterRange range(resultGrid->tree().cbeginLeaf());
    tbb::parallel_for(range, proc);

    return resultGrid;
}





struct GradientXComputer {
    using GridType = openvdb::FloatGrid;
    using TreeType = GridType::TreeType;
    using LeafNode = TreeType::LeafNodeType;
    using IterRange = openvdb::tree::IteratorRange<typename TreeType::LeafCIter>;

    GridType::ConstPtr inputGrid;  // Using ConstPtr for--- input grid
    tbb::enumerable_thread_specific<GridType::Ptr>& threadLocalGrids;  // Thread-local grids

    //constructor
    GradientXComputer(GridType::ConstPtr ig, tbb::enumerable_thread_specific<GridType::Ptr>& tlg)
        : inputGrid(ig), threadLocalGrids(tlg) {}

    void operator()(IterRange& range) const {
        auto& gradXGrid = threadLocalGrids.local();
        GridType::Accessor accessor = gradXGrid->getAccessor();  // local accessor __for each thread__ to write data

        for (; range; ++range) {
            const LeafNode& leaf = *range.iterator();
            for (auto iter = leaf.beginValueOn(); iter; ++iter) {
                const openvdb::Coord xyz = iter.getCoord();
                float valueLeft = inputGrid->tree().isValueOn(xyz.offsetBy(-1, 0, 0)) ?
                    inputGrid->tree().getValue(xyz.offsetBy(-1, 0, 0)) : iter.getValue();
                float valueRight = inputGrid->tree().isValueOn(xyz.offsetBy(1, 0, 0)) ?
                    inputGrid->tree().getValue(xyz.offsetBy(1, 0, 0)) : iter.getValue();
                float gradX = (valueRight - valueLeft) / 2.0f;

                accessor.setValueOn(xyz, gradX);  // using local accessor--write 
            }
        }
    }
};

openvdb::FloatGrid::Ptr calculateGradientX_tbb(const openvdb::FloatGrid::Ptr& grid) {
    openvdb::FloatGrid::Ptr gradXGrid = openvdb::FloatGrid::create();
    gradXGrid->setTransform(grid->transform().copy());

    // create thread-local grids
    tbb::enumerable_thread_specific<openvdb::FloatGrid::Ptr> threadLocalGrids([&]() {
        auto localGrid = openvdb::FloatGrid::create();
        localGrid->setTransform(grid->transform().copy());
        return localGrid;
        });

    GradientXComputer::IterRange range(grid->tree().cbeginLeaf());
    GradientXComputer proc(grid, threadLocalGrids);

    tbb::parallel_for(range, proc);

    // merge the thread-local grids into the final output grid
    for (auto& localGrid : threadLocalGrids) {
        gradXGrid->tree().merge(localGrid->tree());
    }

    return gradXGrid;
}





struct GradientYProcessor {
    using GridType = openvdb::FloatGrid;
    using TreeType = GridType::TreeType;
    using LeafNode = TreeType::LeafNodeType;
    using IterRange = openvdb::tree::IteratorRange<typename TreeType::LeafCIter>;

    GridType::ConstPtr inputGrid;  // Using ConstPtr for input grid [using ConstPtr for the input grid to ensure it remains read-only]
    tbb::enumerable_thread_specific<GridType::Ptr>& threadLocalGrids;  // Thread-local grids

    GradientYProcessor(GridType::ConstPtr ig, tbb::enumerable_thread_specific<GridType::Ptr>& tlg)
        : inputGrid(ig), threadLocalGrids(tlg) {}

    void operator()(IterRange& range) const {
        auto& gradYGrid = threadLocalGrids.local();
        GridType::Accessor accessor = gradYGrid->getAccessor();  // Local accessor for each thread to write data

        for (; range; ++range) {
            const LeafNode& leaf = *range.iterator();
            for (auto iter = leaf.beginValueOn(); iter; ++iter) {
                const openvdb::Coord xyz = iter.getCoord();
                float valueDown = inputGrid->tree().isValueOn(xyz.offsetBy(0, -1, 0)) ?
                    inputGrid->tree().getValue(xyz.offsetBy(0, -1, 0)) : iter.getValue();
                float valueUp = inputGrid->tree().isValueOn(xyz.offsetBy(0, 1, 0)) ?
                    inputGrid->tree().getValue(xyz.offsetBy(0, 1, 0)) : iter.getValue();
                float gradY = (valueUp - valueDown) / 2.0f;

                accessor.setValueOn(xyz, gradY);  // Write using local accessor
            }
        }
    }
};

openvdb::FloatGrid::Ptr calculateGradientY_tbb(const openvdb::FloatGrid::Ptr& grid) {
    openvdb::FloatGrid::Ptr gradYGrid = openvdb::FloatGrid::create();
    gradYGrid->setTransform(grid->transform().copy());

    // Create thread-local grids
    tbb::enumerable_thread_specific<openvdb::FloatGrid::Ptr> threadLocalGrids([&]() {
        auto localGrid = openvdb::FloatGrid::create();
        localGrid->setTransform(grid->transform().copy());
        return localGrid;
        });

    GradientYProcessor::IterRange range(grid->tree().cbeginLeaf());
    GradientYProcessor proc(grid, threadLocalGrids);

    tbb::parallel_for(range, proc);

    // Merge thread-local grids ------- final output grid
    for (auto& localGrid : threadLocalGrids) {
        gradYGrid->tree().merge(localGrid->tree());
    }

    return gradYGrid;
}

struct GradientZProcessor {
    using GridType = openvdb::FloatGrid;
    using TreeType = GridType::TreeType;
    using LeafNode = TreeType::LeafNodeType;
    using IterRange = openvdb::tree::IteratorRange<typename TreeType::LeafCIter>;

    GridType::ConstPtr inputGrid;  // Using ConstPtr for input grid
    tbb::enumerable_thread_specific<GridType::Ptr>& threadLocalGrids;  // Thread-local grids

    GradientZProcessor(GridType::ConstPtr ig, tbb::enumerable_thread_specific<GridType::Ptr>& tlg)
        : inputGrid(ig), threadLocalGrids(tlg) {}

    void operator()(IterRange& range) const {
        auto& gradZGrid = threadLocalGrids.local();
        GridType::Accessor accessor = gradZGrid->getAccessor();  // Local accessor for each thread to write data

        for (; range; ++range) {
            const LeafNode& leaf = *range.iterator();
            for (auto iter = leaf.beginValueOn(); iter; ++iter) {
                const openvdb::Coord xyz = iter.getCoord();
                float valueBelow = inputGrid->tree().isValueOn(xyz.offsetBy(0, 0, -1)) ?
                    inputGrid->tree().getValue(xyz.offsetBy(0, 0, -1)) : iter.getValue();
                float valueAbove = inputGrid->tree().isValueOn(xyz.offsetBy(0, 0, 1)) ?
                    inputGrid->tree().getValue(xyz.offsetBy(0, 0, 1)) : iter.getValue();
                float gradZ = (valueAbove - valueBelow) / 2.0f;

                accessor.setValueOn(xyz, gradZ);  // Write using local accessor
            }
        }
    }
};

openvdb::FloatGrid::Ptr calculateGradientZ_tbb(const openvdb::FloatGrid::Ptr& grid) {
    openvdb::FloatGrid::Ptr gradZGrid = openvdb::FloatGrid::create();
    gradZGrid->setTransform(grid->transform().copy());

    // Create thread-local grids
    tbb::enumerable_thread_specific<openvdb::FloatGrid::Ptr> threadLocalGrids([&]() {
        auto localGrid = openvdb::FloatGrid::create();
        localGrid->setTransform(grid->transform().copy());
        return localGrid;
        });

    GradientZProcessor::IterRange range(grid->tree().cbeginLeaf());
    GradientZProcessor proc(grid, threadLocalGrids);

    tbb::parallel_for(range, proc);

    // Merge
    for (auto& localGrid : threadLocalGrids) {
        gradZGrid->tree().merge(localGrid->tree());
    }

    return gradZGrid;
}


struct GradientMagnitudeComputer {
    using GridType = openvdb::FloatGrid;
    using TreeType = GridType::TreeType;
    using LeafNode = TreeType::LeafNodeType;
    using IterRange = openvdb::tree::IteratorRange<typename TreeType::LeafCIter>;

    GridType::ConstPtr gradXGrid;
    GridType::ConstPtr gradYGrid;
    GridType::ConstPtr gradZGrid;
    tbb::enumerable_thread_specific<GridType::Ptr>& threadLocalGrids;

    GradientMagnitudeComputer(GridType::ConstPtr gx, GridType::ConstPtr gy, GridType::ConstPtr gz, tbb::enumerable_thread_specific<GridType::Ptr>& tlg)
        : gradXGrid(gx), gradYGrid(gy), gradZGrid(gz), threadLocalGrids(tlg) {}

    void operator()(IterRange& range) const {
        auto& gradMagGrid = threadLocalGrids.local();
        GridType::Accessor accessor = gradMagGrid->getAccessor();  // Local accessor

        auto accessorX = gradXGrid->getConstAccessor();
        auto accessorY = gradYGrid->getConstAccessor();
        auto accessorZ = gradZGrid->getConstAccessor();

        for (; range; ++range) {
            const LeafNode& leaf = *range.iterator();
            for (auto iter = leaf.beginValueOn(); iter; ++iter) {
                const openvdb::Coord& xyz = iter.getCoord();
                float gradX = accessorX.getValue(xyz);
                float gradY = accessorY.getValue(xyz);
                float gradZ = accessorZ.getValue(xyz);
                float magnitude = openvdb::math::Sqrt(gradX * gradX + gradY * gradY + gradZ * gradZ);

                accessor.setValueOn(xyz, magnitude);  
            }
        }
    }
};

openvdb::FloatGrid::Ptr calculateGradientMagnitude_tbb(
    const openvdb::FloatGrid::Ptr& gradXGrid,
    const openvdb::FloatGrid::Ptr& gradYGrid,
    const openvdb::FloatGrid::Ptr& gradZGrid) {

    openvdb::FloatGrid::Ptr gradMagGrid = openvdb::FloatGrid::create();
    gradMagGrid->setTransform(gradXGrid->transform().copy());
    gradMagGrid->setGridClass(openvdb::GRID_LEVEL_SET);

    // Create thread-local grids
    tbb::enumerable_thread_specific<openvdb::FloatGrid::Ptr> threadLocalGrids([&]() {
        auto localGrid = openvdb::FloatGrid::create();
        localGrid->setTransform(gradXGrid->transform().copy());
        return localGrid;
        });

    GradientMagnitudeComputer::IterRange range(gradXGrid->tree().cbeginLeaf());
    GradientMagnitudeComputer proc(gradXGrid, gradYGrid, gradZGrid, threadLocalGrids);

    tbb::parallel_for(range, proc);

    // Merge 
    for (auto& localGrid : threadLocalGrids) {
        gradMagGrid->tree().merge(localGrid->tree());
    }

    return gradMagGrid;
}



openvdb::FloatGrid::Ptr dist3D(openvdb::FloatGrid::Ptr phi, openvdb::FloatGrid::Ptr binary_boundary) {
    // Initialize grid dimensions
    const openvdb::CoordBBox bbox = phi->evalActiveVoxelBoundingBox();
    int width = bbox.dim().x();
    int height = bbox.dim().y();
    int length = bbox.dim().z();

    // Create binary boundary grid
    binary_boundary->setTransform(phi->transform().copy());

    // Neighbor offsets
    std::vector<openvdb::Coord> neighbors = {
        openvdb::Coord(0, 0, 1),
        openvdb::Coord(0, 1, 0),
        openvdb::Coord(1, 0, 0),
        openvdb::Coord(-1, 0, 0),
        openvdb::Coord(0, -1, 0),
        openvdb::Coord(0, 0, -1)
    };

    openvdb::FloatGrid::Accessor phiAccessor = phi->getAccessor();
    openvdb::FloatGrid::Accessor binaryAccessor = binary_boundary->getAccessor();

    // Identifying boundary cells using tiles
    for (auto leafIter = phi->tree().cbeginLeaf(); leafIter; ++leafIter) {
        for (auto iter = leafIter->cbeginValueOn(); iter; ++iter) {
            const openvdb::Coord& xyz = iter.getCoord();
            for (const auto& neighbor : neighbors) {
                openvdb::Coord nxyz = xyz + neighbor;
                if (phiAccessor.isValueOn(nxyz) && (iter.getValue() * phiAccessor.getValue(nxyz) <= 0)) {
                    binaryAccessor.setValueOn(xyz);
                    binaryAccessor.setValueOn(nxyz);
                }
            }
        }
    }

    // Create distance grid and initialize with large values
    openvdb::FloatGrid::Ptr dist = openvdb::FloatGrid::create(35.0f);
    dist->setTransform(phi->transform().copy());

    openvdb::FloatGrid::Accessor distAccessor = dist->getAccessor();

    // Calculate distance for all boundary cells to the contour using tiles
    for (auto leafIter = binary_boundary->tree().cbeginLeaf(); leafIter; ++leafIter) {
        for (auto iter = leafIter->cbeginValueOn(); iter; ++iter) {
            const openvdb::Coord& xyz = iter.getCoord();
            for (const auto& neighbor : neighbors) {
                openvdb::Coord nxyz = neighbor + xyz;
                if (binaryAccessor.isValueOn(nxyz)) {
                    float da = std::abs(phiAccessor.getValue(xyz)) / std::abs(phiAccessor.getValue(nxyz) - phiAccessor.getValue(xyz));
                    float db = std::abs(phiAccessor.getValue(nxyz)) / std::abs(phiAccessor.getValue(nxyz) - phiAccessor.getValue(xyz));
                    distAccessor.setValue(xyz, std::min(distAccessor.getValue(xyz), da));
                    distAccessor.setValue(nxyz, std::min(distAccessor.getValue(nxyz), db));
                }
            }
        }
    }

    // Fast sweeping algorithm
    const int NSweeps = 4;
    const int dirX[NSweeps][3] = { {0, width - 1, 1}, {width - 1, 0, -1}, {width - 1, 0, -1}, {0, width - 1, 1} };
    const int dirY[NSweeps][3] = { {0, height - 1, 1}, {0, height - 1, 1}, {height - 1, 0, -1}, {height - 1, 0, -1} };
    const int dirZ[NSweeps][3] = { {0, length - 1, 1}, {0, length - 1, 1}, {0, length - 1, 1}, {0, length - 1, 1} };
    const double h = 1.0, f = 1.0, eps = 1e-6;

    for (int s = 0; s < NSweeps; s++) {
        for (int iy = dirY[s][0]; dirY[s][2] * iy <= dirY[s][1]; iy += dirY[s][2]) {
            for (int ix = dirX[s][0]; dirX[s][2] * ix <= dirX[s][1]; ix += dirX[s][2]) {
                for (int iz = dirZ[s][0]; dirZ[s][2] * iz <= dirZ[s][1]; iz += dirZ[s][2]) {
                    openvdb::Coord xyz(ix, iy, iz);
                    if (phi->tree().isValueOn(xyz)) {
                        if (!binaryAccessor.isValueOn(xyz)) {
                            double aa[3], tmp;
                            if (iy == 0) aa[1] = distAccessor.getValue(openvdb::Coord(ix, iy + 1, iz));
                            else if (iy == height - 1) aa[1] = distAccessor.getValue(openvdb::Coord(ix, iy - 1, iz));
                            else aa[1] = std::min(distAccessor.getValue(openvdb::Coord(ix, iy - 1, iz)),
                                distAccessor.getValue(openvdb::Coord(ix, iy + 1, iz)));

                            if (ix == 0) aa[0] = distAccessor.getValue(openvdb::Coord(ix + 1, iy, iz));
                            else if (ix == width - 1) aa[0] = distAccessor.getValue(openvdb::Coord(ix - 1, iy, iz));
                            else aa[0] = std::min(distAccessor.getValue(openvdb::Coord(ix - 1, iy, iz)),
                                distAccessor.getValue(openvdb::Coord(ix + 1, iy, iz)));

                            if (iz == 0) aa[2] = distAccessor.getValue(openvdb::Coord(ix, iy, iz + 1));
                            else if (iz == length - 1) aa[2] = distAccessor.getValue(openvdb::Coord(ix, iy, iz - 1));
                            else aa[2] = std::min(distAccessor.getValue(openvdb::Coord(ix, iy, iz - 1)),
                                distAccessor.getValue(openvdb::Coord(ix, iy, iz + 1)));

                            // Simple bubble sort
                            if (aa[0] > aa[1]) { tmp = aa[0]; aa[0] = aa[1]; aa[1] = tmp; }
                            if (aa[1] > aa[2]) { tmp = aa[1]; aa[1] = aa[2]; aa[2] = tmp; }
                            if (aa[0] > aa[1]) { tmp = aa[0]; aa[0] = aa[1]; aa[1] = tmp; }

                            double d_curr = aa[0] + h * f;
                            double d_new;
                            if (d_curr <= (aa[1] + eps)) {
                                d_new = d_curr;
                            }
                            else {
                                double a = 2.0, b = -2.0 * (aa[0] + aa[1]), c = aa[0] * aa[0] + aa[1] * aa[1] - h * h * f * f;
                                double D = sqrt(b * b - 4.0 * a * c);
                                d_curr = ((-b + D) > (-b - D) ? (-b + D) : (-b - D)) / (2.0 * a);

                                if (d_curr <= (aa[2] + eps)) {
                                    d_new = d_curr;
                                }
                                else {
                                    a = 3.0;
                                    b = -2.0 * (aa[0] + aa[1] + aa[2]);
                                    c = aa[0] * aa[0] + aa[1] * aa[1] + aa[2] * aa[2] - h * h * f * f;
                                    D = sqrt(b * b - 4.0 * a * c);
                                    d_new = ((-b + D) > (-b - D) ? (-b + D) : (-b - D)) / (2.0 * a);
                                }
                            }
                            distAccessor.setValue(xyz, std::min(distAccessor.getValue(xyz), float(d_new)));
                        }
                    }
                }
            }
        }
    }

    return dist;
}



openvdb::FloatGrid::Ptr dist2sdf3D(openvdb::FloatGrid::Ptr distance, openvdb::FloatGrid::Ptr binary_boundary) {
    // Initialize grid dimensions
    const openvdb::CoordBBox bbox = distance->evalActiveVoxelBoundingBox();
    int width = bbox.dim().x();
    int height = bbox.dim().y();
    int length = bbox.dim().z();

    std::cout << "start_sdf" << "\n";
    // Initialize SDF and frozen cells grids
    auto SDF = openvdb::FloatGrid::create(std::numeric_limits<float>::max());
    auto frozenCells = openvdb::Int32Grid::create(0);

    SDF->setTransform(distance->transform().copy());
    frozenCells->setTransform(distance->transform().copy());

    // Copy distance values to SDF grid
    for (auto iter = distance->cbeginValueOn(); iter; ++iter) {
        const openvdb::Coord& xyz = iter.getCoord();
        SDF->tree().setValue(xyz, iter.getValue());
    }
    std::cout << "start_sdf(2)" << "\n";
    // Fill the frozenCells grid for active voxels in binary_boundary
    for (auto iter = binary_boundary->cbeginValueOn(); iter; ++iter) {
        const openvdb::Coord& xyz = iter.getCoord();
        if (distance->tree().isValueOn(xyz)) {  // Check if the voxel is active in the distance grid
            frozenCells->tree().setValue(xyz, 1);
        }
    }

    // Turn the whole input distance field to negative
    for (auto iter = distance->cbeginValueOn(); iter; ++iter) {
        const openvdb::Coord& xyz = iter.getCoord();
        float value = iter.getValue();
        SDF->tree().setValue(xyz, -value);
    }
    std::cout << "start_sdf(3)" << "\n";
    std::queue<openvdb::Coord> queue;

    // Initialize the queue with all active voxels in distance that are not frozen
    for (auto iter = distance->cbeginValueOn(); iter; ++iter) {
        const openvdb::Coord& xyz = iter.getCoord();
        if (frozenCells->tree().getValue(xyz) == 0) {
            queue.push(xyz);
        }
    }

    std::cout << "start_sdf(4)" << "\n";

    // Define a list of neighbor coordinate offsets
    const std::vector<openvdb::Coord> offsets = {
        openvdb::Coord(-1, 0, 0), openvdb::Coord(1, 0, 0),
        openvdb::Coord(0, -1, 0), openvdb::Coord(0, 1, 0),
        openvdb::Coord(0, 0, -1), openvdb::Coord(0, 0, 1)
    };

    //// A simple voxel flood using BFS
    /*while (!queue.empty()) {
        openvdb::Coord xyz = queue.front();
        queue.pop();

        if (frozenCells->tree().getValue(xyz) == 1) continue;

        float val = -1.0f * SDF->tree().getValue(xyz);
        SDF->tree().setValue(xyz, val);
        frozenCells->tree().setValue(xyz, 1);


        for (const auto& offset : offsets) {
            openvdb::Coord neighbor = xyz + offset;
            if (distance->tree().isValueOn(neighbor) && frozenCells->tree().getValue(neighbor) == 0) {
                queue.push(neighbor);
            }
        }
    }*/

    const int VOXEL_FROZEN = 1;
    const int VOXEL_IN_QUEUE = 2;

    // A simple voxel flood using BFS
    while (!queue.empty()) {
        openvdb::Coord xyz = queue.front();
        queue.pop();
        if (frozenCells->tree().getValue(xyz) != VOXEL_FROZEN) {
            float val = -1.0f * SDF->tree().getValue(xyz);
            SDF->tree().setValue(xyz, val);
            frozenCells->tree().setValue(xyz, VOXEL_FROZEN); // Freeze cell when done

            // Check all neighbors
            for (const auto& offset : offsets) {
                openvdb::Coord neighbor = xyz + offset;
                if (distance->tree().isValueOn(neighbor)) {
                    int neighborState = frozenCells->tree().getValue(neighbor);
                    if (neighborState == 0) { // Not processed or in queue
                        frozenCells->tree().setValue(neighbor, VOXEL_IN_QUEUE);
                        queue.push(neighbor);
                    }
                }
            }
        }
    }

    std::cout << "start_sdf(6)" << "\n";
    // Update the distance grid with the new SDF values
    for (auto iter = distance->beginValueOn(); iter; ++iter) {
        const openvdb::Coord& xyz = iter.getCoord();
        float value = SDF->tree().getValue(xyz);
        distance->tree().setValue(xyz, value);
    }

    return distance;
}



template<class GridType>
void printNaNValues(GridType& grid) {
    using ValueT = typename GridType::ValueType;

    // Iterate through all active voxels in the grid
    for (auto iter = grid.beginValueOn(); iter; ++iter) {
        const ValueT value = iter.getValue();
        if (std::isnan(value)) {
            openvdb::Coord coord = iter.getCoord();
            std::cout << "NaN value found at (" << coord.x() << ", " << coord.y() << ", " << coord.z() << ")\n";
        }
    }
}



int main()
{
    openvdb::initialize();

    // Load sdf grid from file
    openvdb::io::File file_sdf("C:/openvdb_drop/bin/center_thresh.vdb");
    file_sdf.open();
    openvdb::GridBase::Ptr baseGrid_sdf;
    for (openvdb::io::File::NameIterator nameIter = file_sdf.beginName(); nameIter != file_sdf.endName(); ++nameIter) {
        if (nameIter.gridName() == "LevelSetSphere") {
            baseGrid_sdf = file_sdf.readGrid(nameIter.gridName());
            break;
        }
    }
    file_sdf.close();

    // Load input grid from file
    openvdb::io::File file_input("C:/openvdb_drop/bin/bmp_thresh.vdb");
    file_input.open();
    openvdb::GridBase::Ptr baseGrid_input;
    for (openvdb::io::File::NameIterator nameIter = file_input.beginName(); nameIter != file_input.endName(); ++nameIter) {
        if (nameIter.gridName() == "LevelSetSphere") {
            baseGrid_input = file_input.readGrid(nameIter.gridName());
            break;
        }
    }
    file_input.close();


    // grids to FloatGrid
    openvdb::FloatGrid::Ptr sdf_grid = openvdb::gridPtrCast<openvdb::FloatGrid>(baseGrid_sdf);
    openvdb::FloatGrid::Ptr input_grid = openvdb::gridPtrCast<openvdb::FloatGrid>(baseGrid_input);

    /*openvdb::FloatGrid::Ptr binary_boundary = openvdb::FloatGrid::create();

    sdf_grid1 = dist3D_tree(sdf_grid1, binary_boundary);
    sdf_grid1 = dist2sdf3D(sdf_grid1, binary_boundary);*/

    size_t memUsage = sdf_grid->memUsage();
    std::cout << "Memory usage: " << memUsage << " bytes" << std::endl;

    //exit(1);

    ///////////////// initializing with empty grid/////////////////////

    /*
    
    // grids to FloatGrid
    openvdb::FloatGrid::Ptr sdf_grid1 = openvdb::gridPtrCast<openvdb::FloatGrid>(baseGrid_sdf);
    openvdb::FloatGrid::Ptr input_grid1 = openvdb::gridPtrCast<openvdb::FloatGrid>(baseGrid_input);
    
    float threshold = 110.0f;
   

    openvdb::FloatGrid::Ptr empty_grid = openvdb::FloatGrid::create();
    openvdb::FloatGrid::Ptr sdf_grid = openvdb::FloatGrid::create(1);
    openvdb::FloatGrid::Ptr input_grid = openvdb::FloatGrid::create(140);
    openvdb::FloatGrid::Ptr sdf_grid2 = openvdb::FloatGrid::create();
    
    for (openvdb::FloatGrid::ValueOnIter iter = input_grid1->beginValueOn(); iter; ++iter) {
        if (iter.getValue() <= threshold) {
            openvdb::Coord ijk = iter.getCoord();


            empty_grid->tree().setValueOn(ijk);
        }
    }

   
    openvdb::tools::dilateActiveValues(empty_grid->tree(), 3);

    for (openvdb::FloatGrid::ValueOnCIter iter = input_grid1->cbeginValueOn(); iter; ++iter) {
        openvdb::Coord ijk = iter.getCoord();

        if (empty_grid->tree().isValueOn(ijk)) {

            float value = input_grid1->tree().getValue(ijk);
            input_grid->tree().setValue(ijk, value);
        }

    }

    for (openvdb::FloatGrid::ValueOnCIter iter = input_grid->cbeginValueOn(); iter; ++iter) {
        openvdb::Coord ijk = iter.getCoord();

        sdf_grid->tree().setValueOn(ijk);
    }


    for (openvdb::FloatGrid::ValueOnCIter iter = sdf_grid->cbeginValueOn(); iter; ++iter) {
        openvdb::Coord ijk = iter.getCoord();

        if (sdf_grid1->tree().isValueOn(ijk)) {

            float value = sdf_grid1->tree().getValue(ijk);
            sdf_grid->tree().setValue(ijk, value);
        }

    }

   ///FOT HD IN

    for(openvdb::FloatGrid::ValueOnCIter iter = sdf_grid->cbeginValueOn(); iter; ++iter) {
        openvdb::Coord ijk = iter.getCoord();

        if (sdf_grid->tree().isValueOn(ijk)) {

            float value = sdf_grid->tree().getValue(ijk);
            sdf_grid2->tree().setValue(ijk, value);
        }

    }


    openvdb::FloatGrid::Ptr binary_boundary = openvdb::FloatGrid::create(0.0f);



    sdf_grid = dist3D(sdf_grid, binary_boundary);

    for (openvdb::FloatGrid::ValueOnCIter iter = sdf_grid->cbeginValueOn(); iter; ++iter) {
        openvdb::Coord ijk = iter.getCoord();

        if (empty_grid->tree().isValueOn(ijk)) {

            float value = sdf_grid->tree().getValue(ijk);
            sdf_grid->tree().setValue(ijk, value);
        }

    }

    for (openvdb::FloatGrid::ValueOnCIter iter = binary_boundary->cbeginValueOn(); iter; ++iter) {
        openvdb::Coord ijk = iter.getCoord();

        if (empty_grid->tree().isValueOn(ijk)) {

            float value = binary_boundary->tree().getValue(ijk);
            binary_boundary->tree().setValue(ijk, value);
        }

    }

    */


    openvdb::tools::changeBackground(sdf_grid->tree(), 255);

    for (openvdb::FloatGrid::ValueOnCIter iter = input_grid->cbeginValueOn(); iter; ++iter) {
        sdf_grid->tree().setValueOn(iter.getCoord());
    }


    openvdb::FloatGrid::Ptr binary_boundary = openvdb::FloatGrid::create(0.0f);


    openvdb::tools::changeBackground(sdf_grid->tree(), 0);


    sdf_grid = dist3D(sdf_grid, binary_boundary);

    std::cout << "done" << "\n";

    sdf_grid = dist2sdf3D(sdf_grid, binary_boundary);

    std::cout << "done sdf" << "\n";
    

    openvdb::tools::changeBackground(input_grid->tree(), 255);
    openvdb::tools::changeBackground(sdf_grid->tree(), 1);
    openvdb::FloatGrid::Ptr sdf_grid2 = sdf_grid->deepCopy();
    openvdb::tools::changeBackground(sdf_grid2->tree(), 0);

    input_grid->tree().prune();
    sdf_grid->tree().prune();

    
    
    std::cout << "Total number of sdf voxel: " << sdf_grid->activeVoxelCount() << std::endl;
    //std::cout << "Total number of input voxel: " << input_grid->activeVoxelCount() << std::endl;


    //openvdb::Coord testCoord(60, 60, 60); 
    //float backgroundvalue = 140.0;
    float backgroundvalue = 255.0;


    //float value = accessor.getValue(testCoord);
    //std::cout << "Value at " << testCoord << ": " << value << std::endl;
    openvdb::FloatGrid::Ptr Eout = openvdb::FloatGrid::create(backgroundvalue);
    openvdb::FloatGrid::Ptr Ein = openvdb::FloatGrid::create(backgroundvalue);
    openvdb::FloatGrid::Ptr DIVGrid = openvdb::FloatGrid::create(backgroundvalue);

    openvdb::FloatGrid::Ptr d2PHI_dx2_1 = openvdb::FloatGrid::create(backgroundvalue);
    openvdb::FloatGrid::Ptr d2PHI_dy2_1 = openvdb::FloatGrid::create(backgroundvalue);
    openvdb::FloatGrid::Ptr d2PHI_dz2_1 = openvdb::FloatGrid::create(backgroundvalue);

    openvdb::FloatGrid::Ptr norm_DIVgrid = openvdb::FloatGrid::create(backgroundvalue);


    int T = 13;
    for (int t = 0; t < T; t++) {


        // Start timing
        auto start = std::chrono::high_resolution_clock::now();

        ////Apply the heaviside function
        openvdb::FloatGrid::Ptr HD_out = sdf_grid->deepCopy();
        heaviside_tbb(HD_out);

        
        
        ////Apply the derivetive of heaviside function
        openvdb::FloatGrid::Ptr deri_heavi = sdf_grid->deepCopy();
        deri_heaviside_tbb(deri_heavi);

        

        // multiplying  heaviside_image with input_image
       

        openvdb::FloatGrid::Ptr I_out = multiplyGrids_tbb(input_grid, HD_out, backgroundvalue);

        


        openvdb::FloatGrid::Ptr HD_in = sdf_grid2->deepCopy();
        heaviside_tbb(HD_in);

        
        // inverse the vdb file
        vdb_inverse_tbb(*HD_in);

        // multiplying  inverse_heaviside_image with input_image

        float b_hd_in = 0.0f;
        openvdb::FloatGrid::Ptr I_in = multiplyGrids_tbb(input_grid, HD_in, b_hd_in);
        openvdb::tools::changeBackground(I_in->tree(),b_hd_in);


        float sigma = 3.0f;

        // Apply Gaussian filter
        openvdb::FloatGrid::Ptr I_out_blurred = I_out->deepCopy();
        openvdb::FloatGrid::Ptr HD_out_blurred = HD_out->deepCopy();
        openvdb::tools::Filter<openvdb::FloatGrid> filter1(*I_out_blurred);
        filter1.setGrainSize(3);
        openvdb::tools::Filter<openvdb::FloatGrid> filter3(*HD_out_blurred);
        filter3.setGrainSize(3);
        
       

        
        filter1.gaussian(sigma);
        filter3.gaussian(sigma);


        // Apply Gaussian filter
        openvdb::FloatGrid::Ptr I_in_blurred = I_in->deepCopy();
        openvdb::FloatGrid::Ptr HD_in_blurred = HD_in->deepCopy();
        openvdb::tools::Filter<openvdb::FloatGrid> filter2(*I_in_blurred);
        filter2.setGrainSize(3);
        openvdb::tools::Filter<openvdb::FloatGrid> filter4(*HD_in_blurred);
        filter4.setGrainSize(3);

        filter2.gaussian(sigma);
        filter4.gaussian(sigma);

        

        //fout
        openvdb::FloatGrid::Ptr f_out = divideGrids_tbb(I_out_blurred, HD_out_blurred);

        

        //fin
        openvdb::FloatGrid::Ptr f_in = divideGrids_tbb(I_in_blurred, HD_in_blurred );

        openvdb::tools::changeBackground(f_in->tree(), backgroundvalue);
        
    
       
         //DIVGrid and normDiv_Grid
        
        
        //  x and y and z of the gradient

        float backvalueGrad = 1.0f;
        
        openvdb::FloatGrid::Ptr dPHI_dx = calculateGradientX_tbb(sdf_grid);
        //openvdb::tools::changeBackground(dPHI_dx->tree(), backvalueGrad);
        openvdb::FloatGrid::Ptr dPHI_dy = calculateGradientY_tbb(sdf_grid);
        //openvdb::tools::changeBackground(dPHI_dy->tree(), backvalueGrad);
        openvdb::FloatGrid::Ptr dPHI_dz = calculateGradientZ_tbb(sdf_grid);
        //openvdb::tools::changeBackground(dPHI_dz->tree(), backvalueGrad);

       
        // Calculate the gradient magnitude grid
        openvdb::FloatGrid::Ptr gradMag = calculateGradientMagnitude_tbb(dPHI_dx, dPHI_dy, dPHI_dz);
        //openvdb::tools::changeBackground(gradMag->tree(), backvalueGrad);

       

        //  x and y and z of the gradient
        openvdb::FloatGrid::Ptr d2PHI_dx_n_n = calculateGradientX_tbb(dPHI_dx);
        //openvdb::tools::changeBackground(d2PHI_dx_n_n->tree(), backvalueGrad);
        openvdb::FloatGrid::Ptr d2PHI_dy_n_n = calculateGradientY_tbb(dPHI_dy);
        //openvdb::tools::changeBackground(d2PHI_dy_n_n->tree(), backvalueGrad);
        openvdb::FloatGrid::Ptr d2PHI_dz_n_n = calculateGradientZ_tbb(dPHI_dz);
        //openvdb::tools::changeBackground(d2PHI_dz_n_n->tree(), backvalueGrad);

        
        openvdb::FloatGrid::Accessor DIVAccessor = DIVGrid->getAccessor();
        openvdb::FloatGrid::Accessor d2PHI_dx_n_nAccessor = d2PHI_dx_n_n->getAccessor();
        openvdb::FloatGrid::Accessor d2PHI_dy_n_nAccessor = d2PHI_dy_n_n->getAccessor();
        openvdb::FloatGrid::Accessor d2PHI_dz_n_nAccessor = d2PHI_dz_n_n->getAccessor();

        // active values of the first grid
        for (openvdb::FloatGrid::ValueOnIter iter = d2PHI_dy_n_n->beginValueOn(); iter; ++iter) {
            const openvdb::Coord& coord = iter.getCoord();

            // Multiply 
            DIVAccessor.setValue(coord, d2PHI_dx_n_nAccessor.getValue(coord) + d2PHI_dy_n_nAccessor.getValue(coord) + d2PHI_dz_n_nAccessor.getValue(coord));
        }

       

        d2PHI_dx2_1 = divideGrids_tbb(dPHI_dx, gradMag);
        d2PHI_dy2_1 = divideGrids_tbb(dPHI_dy, gradMag);
        d2PHI_dz2_1 = divideGrids_tbb(dPHI_dz, gradMag);

        openvdb::FloatGrid::Ptr d2PHI_dx2 = calculateGradientX_tbb(d2PHI_dx2_1);
        //openvdb::tools::changeBackground(d2PHI_dx2->tree(), backvalueGrad);
        openvdb::FloatGrid::Ptr d2PHI_dy2 = calculateGradientY_tbb(d2PHI_dy2_1);
        //openvdb::tools::changeBackground(d2PHI_dy2->tree(), backvalueGrad);
        openvdb::FloatGrid::Ptr d2PHI_dz2 = calculateGradientZ_tbb(d2PHI_dz2_1);
        //openvdb::tools::changeBackground(d2PHI_dz2->tree(), backvalueGrad);

        openvdb::FloatGrid::Accessor norm_DIVgridAccessor = norm_DIVgrid->getAccessor();
        openvdb::FloatGrid::Accessor d2PHI_dx2Accessor = d2PHI_dx2->getAccessor();
        openvdb::FloatGrid::Accessor d2PHI_dy2Accessor = d2PHI_dy2->getAccessor();
        openvdb::FloatGrid::Accessor d2PHI_dz2Accessor = d2PHI_dz2->getAccessor();

        for (openvdb::FloatGrid::ValueOnIter iter = d2PHI_dy2->beginValueOn(); iter; ++iter) {
            const openvdb::Coord& coord = iter.getCoord();

            // Multiply 
            norm_DIVgridAccessor.setValue(coord, d2PHI_dx2Accessor.getValue(coord) + d2PHI_dy2Accessor.getValue(coord) + d2PHI_dz2Accessor.getValue(coord));
        }

        
        
        //Eout and Ein Calculation
        
        //float sigma = 3.0f; 
        int sigma1 = 0;

       
        
        parallelVoxelComputation(input_grid, f_out, f_in, Eout, Ein, sigma1);

        
       

        float meu = 255 * 255 * 0.0001;
        float dt = 0.1; //0.5
        float f = 0.1; //1.5
        openvdb::FloatGrid::Accessor Eoutaccessor = Eout->getAccessor();
        openvdb::FloatGrid::Accessor Einaccessor = Ein->getAccessor();
        openvdb::FloatGrid::Accessor sdf_accessor = sdf_grid->getAccessor();
        openvdb::FloatGrid::Accessor deri_heavi_accessor = deri_heavi->getAccessor();
        openvdb::FloatGrid::Accessor DivGrid_accessor = DIVGrid->getAccessor();
        openvdb::FloatGrid::Accessor normDivGrid_accessor = norm_DIVgrid->getAccessor();

        float lambda_out = 1.3f;
        

        // openvdb::FloatGrid::Ptr sdf_grid_new = sdf_grid->deepCopy();

        for (openvdb::FloatGrid::ValueOnIter iter = input_grid->beginValueOn(); iter; ++iter) {
            const openvdb::Coord& coord = iter.getCoord();
           // sdf_accessor.setValue(coord, sdf_accessor.getValue(coord) - (dt * 1.3 * (deri_heavi_accessor.getValue(coord) * (Eoutaccessor.getValue(coord) - Einaccessor.getValue(coord)))) );
            sdf_accessor.setValue(coord, sdf_accessor.getValue(coord) - (dt * f * (deri_heavi_accessor.getValue(coord) * (Eoutaccessor.getValue(coord) - Einaccessor.getValue(coord)))) + (dt * meu * deri_heavi_accessor.getValue(coord) * normDivGrid_accessor.getValue(coord)) + (dt * DivGrid_accessor.getValue(coord) - dt * normDivGrid_accessor.getValue(coord)));
        }

    


        auto finish = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = finish - start;
        std::cout << "Elapsed time: " << elapsed.count() << " seconds" << std::endl;

        printNaNValues(*sdf_grid);


    }
    for (openvdb::FloatGrid::ValueOnIter iter = sdf_grid->beginValueOn(); iter; ++iter) {
        if (*iter > 0.0f) {
            //if its value is greater than zero
            iter.setValueOff();
        }
    }
    //exit(1);
    std::string phi = "C:/openvdb_drop/bin/skfm_phi.vdb";
    openvdb::initialize();

    // Create a VDB file object.
    openvdb::io::File fileEout(phi);
    fileEout.write({ sdf_grid });
    //// Close the file. 
    fileEout.close();

    tira::volume<float> I2(200, 200, 200);
    vdb2img3D(*sdf_grid, I2); 
    I2.save_npy("C:/Users/meher/spyder/HD13.npy");

    return 0;
}
