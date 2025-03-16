#include "target_image_widget.h"
//#include "seamless.h"

#include <cmath>

namespace USTC_CG
{
using uchar = unsigned char;

TargetImageWidget::TargetImageWidget(
    const std::string& label,
    const std::string& filename)
    : ImageWidget(label, filename)
{
    if (data_)
        back_up_ = std::make_shared<Image>(*data_);
}

void TargetImageWidget::draw()
{
    // Draw the image
    ImageWidget::draw();
    // Invisible button for interactions
    ImGui::SetCursorScreenPos(position_);
    ImGui::InvisibleButton(
        label_.c_str(),
        ImVec2(
            static_cast<float>(image_width_),
            static_cast<float>(image_height_)),
        ImGuiButtonFlags_MouseButtonLeft);
    bool is_hovered_ = ImGui::IsItemHovered();
    // When the mouse is clicked or moving, we would adapt clone function to
    // copy the selected region to the target.

    if (is_hovered_ && ImGui::IsMouseClicked(ImGuiMouseButton_Left))
    {
        mouse_click_event();
    }
    mouse_move_event();
    if (!ImGui::IsMouseDown(ImGuiMouseButton_Left))
    {
        mouse_release_event();
    }
}

void TargetImageWidget::set_source(std::shared_ptr<SourceImageWidget> source)
{
    source_image_ = source;
}

void TargetImageWidget::set_realtime(bool flag)
{
    flag_realtime_updating = flag;
}

void TargetImageWidget::restore()
{
    *data_ = *back_up_;
    update();
}

void TargetImageWidget::set_paste()
{
    clone_type_ = kPaste;
}

void TargetImageWidget::set_seamless()
{
    clone_type_ = kSeamless;
}
void TargetImageWidget::set_mixing()
{
    clone_type_ = kMixing;
}

void TargetImageWidget::clone()
{
    // The implementation of different types of cloning
    // HW3_TODO: 
    // 1. In this function, you should at least implement the "seamless"
    // cloning labeled by `clone_type_ ==kSeamless`.
    //
    // 2. It is required to improve the efficiency of your seamless cloning to
    // achieve real-time editing. (Use decomposition of sparse matrix before
    // solve the linear system). The real-time updating (update when the mouse
    // is moving) is only available when the checkerboard is selected. 
    if (data_ == nullptr || source_image_ == nullptr ||
        source_image_->get_region_mask() == nullptr)
        return;
    // The selected region in the source image, this would be a binary mask.
    // The **size** of the mask should be the same as the source image.
    // The **value** of the mask should be 0 or 255: 0 for the background and
    // 255 for the selected region.
    std::shared_ptr<Image> mask = source_image_->get_region_mask();
    

    switch (clone_type_)
    {
        case USTC_CG::TargetImageWidget::kDefault: break;
        case USTC_CG::TargetImageWidget::kPaste:
        {
            restore();

            for (int x = 0; x < mask->width(); ++x)
            {
                for (int y = 0; y < mask->height(); ++y)
                {
                    int tar_x =
                        static_cast<int>(mouse_position_.x) + x -
                        static_cast<int>(source_image_->get_position().x);
                    int tar_y =
                        static_cast<int>(mouse_position_.y) + y -
                        static_cast<int>(source_image_->get_position().y);
                    if (0 <= tar_x && tar_x < image_width_ && 0 <= tar_y &&
                        tar_y < image_height_ && mask->get_pixel(x, y)[0] > 0)
                    {
                        data_->set_pixel(
                            tar_x,
                            tar_y,
                            source_image_->get_data()->get_pixel(x, y));
                    }
                }
            }//be careful : the size of the mask should be the same as the source image!!!
            break;
        }
        case USTC_CG::TargetImageWidget::kSeamless:
        {
            restore();
            std::pair<int , int> p(0,0), q(0,0);
            for(int i = 0; i < mask->width(); i++)
            {
                for(int j = 0; j < mask->height(); j++)
                {
                    if(mask->get_pixel(i, j)[0] > 0)
                    {
                        p.first = i;
                        p.second = j;
                        i = mask->width();
                        j = mask->height();
                    }
                }
            
            }
            for(int i = mask->width() - 1; i >= 0; i--)
            {
                for(int j = mask->height() - 1; j >= 0; j--)
                {
                    if(mask->get_pixel(i, j)[0] > 0)
                    {
                        q.first = i;
                        q.second = j;
                        i = -1;
                        j = -1;
                    }
                }
            }   
            
            // seamless_clone_.set_size(q.first-p.first+1, q.second-p.second+1);
            // seamless_clone_.set_offset(
            //     static_cast<int>(mouse_position_.x) - static_cast<int>(source_image_->get_position().x) + p.first,
            //     static_cast<int>(mouse_position_.y) - static_cast<int>(source_image_->get_position().y) + p.second
            // );
            seamless_clone_.set_images(source_image_->get_data(), data_, mask);
            seamless_clone_.set_right_to_mouse(static_cast<int>(mouse_position_.x-source_image_->get_position().x),static_cast<int>(mouse_position_.y-source_image_->get_position().y));
            if(!seamless_clone_.judge_A_precomputed()){
                seamless_clone_.precompute_A();
            }
            //seamless_clone_.precompute_A();
            
            std::shared_ptr<Image> result = seamless_clone_.solve();
            //source_image_->get_data() is the source image, data_ is the target image, mask is the mask of the source image
            //source is not a image class so it must use the get data founction to get the image data.
            //source_image_ is a address
            break;
        }
        case USTC_CG::TargetImageWidget::kMixing:
        {
            restore();
            
            std::pair<int , int> p(0,0), q(0,0);
            for(int i = 0; i < mask->width(); i++)
            {
                for(int j = 0; j < mask->height(); j++)
                {
                    if(mask->get_pixel(i, j)[0] > 0)
                    {
                        p.first = i;
                        p.second = j;
                        i = mask->width();
                        j = mask->height();
                    }
                }
            
            }
            for(int i = mask->width() - 1; i >= 0; i--)
            {
                for(int j = mask->height() - 1; j >= 0; j--)
                {
                    if(mask->get_pixel(i, j)[0] > 0)
                    {
                        q.first = i;
                        q.second = j;
                        i = -1;
                        j = -1;
                    }
                }
            }

            // USTC_CG::SeamlessClone clone;
            mixing_clone_.set_images(source_image_->get_data(), data_, mask);
            // mixing_clone_.set_size(q.first-p.first+1, q.second-p.second+1);
            // mixing_clone_.set_offset(
            //     static_cast<int>(mouse_position_.x) - static_cast<int>(source_image_->get_position().x) + p.first,
            //     static_cast<int>(mouse_position_.y) - static_cast<int>(source_image_->get_position().y) + p.second
            // );
            
            if(!mixing_clone_.judge_A_precomputed()){
                mixing_clone_.precompute_A();
            }
            mixing_clone_.set_right_to_mouse(static_cast<int>(mouse_position_.x) - static_cast<int>(source_image_->get_position().x),static_cast<int>(mouse_position_.y) - static_cast<int>(source_image_->get_position().y));

            //clone.precompute_A();
            mixing_clone_.set_images(source_image_->get_data(), data_, mask);
            std::shared_ptr<Image> result = mixing_clone_.solve();
            //source_image_->get_data() is the source image, data_ is the target image, mask is the mask of the source image
            //source is not a image class so it must use the get data founction to get the image data.
            //source_image_ is a address
            break;
        }
        default: break;
    }

    update();
}

void TargetImageWidget::mouse_click_event()
{
    edit_status_ = true;
    mouse_position_ = mouse_pos_in_canvas();
    clone();
}

void TargetImageWidget::mouse_move_event()
{
    if (edit_status_)
    {
        mouse_position_ = mouse_pos_in_canvas();
        if (flag_realtime_updating)
            clone();
    }
}

void TargetImageWidget::mouse_release_event()
{
    if (edit_status_)
    {
        edit_status_ = false;
    }
}

ImVec2 TargetImageWidget::mouse_pos_in_canvas() const
{
    ImGuiIO& io = ImGui::GetIO();
    return ImVec2(io.MousePos.x - position_.x, io.MousePos.y - position_.y);
}
}  // namespace USTC_CG