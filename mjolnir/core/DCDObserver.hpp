#ifndef MJOLNIR_CORE_DCD_OBSERVER_HPP
#define MJOLNIR_CORE_DCD_OBSERVER_HPP
#include <mjolnir/core/ObserverBase.hpp>
#include <mjolnir/util/progress_bar.hpp>
#include <iostream>
#include <fstream>
#include <iomanip>

namespace mjolnir
{

template<typename traitsT>
class DCDObserver final : public ObserverBase<traitsT>
{
  public:
    using base_type         = ObserverBase<traitsT>;
    using traits_type       = typename base_type::traits_type;
    using real_type         = typename base_type::real_type;
    using coordinate_type   = typename base_type::coordinate_type;
    using system_type       = typename base_type::system_type;
    using forcefield_type   = typename base_type::forcefield_type;
    using progress_bar_type = progress_bar<50>;

  public:

    DCDObserver(const std::string& filename_prefix, bool output_progress = false)
      : base_type(), output_progress_(output_progress), progress_bar_(1),
        prefix_(filename_prefix),
        pos_name_(filename_prefix + std::string("_position.dcd")),
        vel_name_(filename_prefix + std::string("_velocity.dcd")),
        ene_name_ (filename_prefix + std::string(".ene"))
    {
        // clear files and throw an error if the files cannot be opened.
        this->clear_file(this->pos_name_);
        this->clear_file(this->vel_name_);
        this->clear_file(this->ene_name_);
    }
    ~DCDObserver() override = default;

    void initialize(const std::size_t total_step,
                    const system_type& sys, const forcefield_type& ff) override
    {
        this->progress_bar_.reset(total_step); // set total_step

        this->write_header(this->pos_name_, total_step, sys, ff, "CORD");
        this->write_header(this->vel_name_, total_step, sys, ff, "VELO");

        // buffer to convert sys and dcd format
        this->buffer_x_.resize(sys.size());
        this->buffer_y_.resize(sys.size());
        this->buffer_z_.resize(sys.size());

        std::ofstream ofs(this->ene_name_, std::ios::app);
        ofs << "# timestep  " << ff.list_energy_name() << " kinetic_energy\n";
        return;
    }

    void finalize(const std::size_t,
                  const system_type&, const forcefield_type&) override
    {
        // update # of frames in the header region
        {
            std::ofstream ofs(this->pos_name_, std::ios::binary | std::ios::app);
            // skip the first block signature
            ofs.seekp(sizeof(std::int32_t), std::ios::beg);

            const std::int32_t number_of_frames(this->number_of_frames_);
            ofs.write(as_bytes(number_of_frames), sizeof(std::int32_t));
        }
        {
            std::ofstream ofs(this->vel_name_, std::ios::binary | std::ios::app);
            // skip the first block signature
            ofs.seekp(sizeof(std::int32_t), std::ios::beg);

            const std::int32_t number_of_frames(this->number_of_frames_);
            ofs.write(as_bytes(number_of_frames), sizeof(std::int32_t));
        }
        return;
    }

    void output(const std::size_t step,
                const system_type& sys, const forcefield_type& ff) override;

    std::string const& prefix() const noexcept override {return prefix_;}

  private:

    void clear_file(const std::string& fname) const
    {
        std::ofstream ofs(fname);
        if(not ofs.good())
        {
            throw_exception<std::runtime_error>("[error] mjolnir::DCDObserver: "
                    "file open error: ", fname);
        }
        return;
    }

    template<typename T>
    static const char* as_bytes(const T& v) noexcept
    {
        return reinterpret_cast<const char*>(std::addressof(v));
    }

    void write_header(const std::string& fname, const std::size_t total_step_sz,
                      const system_type& sys,   const forcefield_type& ff,
                      const char* signature) const
    {
        std::ofstream ofs(fname, std::ios::binary | std::ios::app);
        if(not ofs.good())
        {
            throw_exception<std::runtime_error>(
                    "[error] mjolnir::DCDObserver: file open error: ", fname);
        }

        /* the first block */
        {
            const std::int32_t block_size(84);
            ofs.write(as_bytes(block_size), sizeof(std::int32_t));
            ofs.write(signature,            4);

            const std::int32_t number_of_frames(0);
            ofs.write(as_bytes(number_of_frames), sizeof(std::int32_t));

            const std::int32_t index_of_first(0);
            ofs.write(as_bytes(index_of_first), sizeof(std::int32_t));

            const std::int32_t save_interval(0);
            ofs.write(as_bytes(save_interval), sizeof(std::int32_t));

            const std::int32_t total_step(total_step_sz);
            ofs.write(as_bytes(total_step), sizeof(std::int32_t));

            const std::int32_t total_chains(1);
            ofs.write(as_bytes(total_chains), sizeof(std::int32_t));

            for(std::size_t i=0; i<4; ++i)
            {
                const std::int32_t zero(0);
                ofs.write(as_bytes(zero), sizeof(std::int32_t));
            }

            const float delta_t(0.0f);
            ofs.write(as_bytes(delta_t), sizeof(float));

            for(std::size_t i=0; i<9; ++i)
            {
                const std::int32_t zero(0);
                ofs.write(as_bytes(zero), sizeof(std::int32_t));
            }

            const std::int32_t version(24);
            ofs.write(as_bytes(version), sizeof(std::int32_t));

            ofs.write(as_bytes(block_size), sizeof(std::int32_t));
        }

        /* the second block */
        {
            const std::int32_t block_size(84);
            ofs.write(as_bytes(block_size), sizeof(std::int32_t));

            const std::int32_t number_of_lines(1);
            ofs.write(as_bytes(number_of_lines), sizeof(std::int32_t));
            const char comment[80] = "Mjolnir -- copyright (c) Toru Niina 2016"
                                     "-now distributed under the MIT License.";
            ofs.write(comment, 80);
            ofs.write(as_bytes(block_size), sizeof(std::int32_t));
        }

        /* the third block */
        {
            const std::int32_t block_size(4);
            ofs.write(as_bytes(block_size), sizeof(std::int32_t));

            const std::int32_t number_of_particles(sys.size());
            ofs.write(as_bytes(number_of_particles), sizeof(std::int32_t));

            ofs.write(as_bytes(block_size), sizeof(std::int32_t));
        }

        return;
    }

    real_type calc_kinetic_energy(const system_type& sys) const
    {
        real_type k = 0.0;
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            k += math::length_sq(sys[i].velocity) * sys[i].mass;
        }
        return k * 0.5;
    }

  private:

    bool output_progress_;
    std::string prefix_;
    std::string pos_name_;
    std::string vel_name_;
    std::string ene_name_;
    std::size_t number_of_frames_;
    std::vector<float> buffer_x_;
    std::vector<float> buffer_y_;
    std::vector<float> buffer_z_;
    progress_bar_type progress_bar_;
};

template<typename traitsT>
inline void DCDObserver<traitsT>::output(
    const std::size_t step, const system_type& sys, const forcefield_type& ff)
{
    number_of_frames_ += 1;
    assert(this->buffer_x_.size() == sys.size());
    assert(this->buffer_y_.size() == sys.size());
    assert(this->buffer_z_.size() == sys.size());
    // ------------------------------------------------------------------------
    // write position
    {
        std::ofstream ofs(this->pos_name_, std::ios::app | std::ios::binary);
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            this->buffer_x_[i] = static_cast<float>(math::X(sys.position(i)));
            this->buffer_y_[i] = static_cast<float>(math::Y(sys.position(i)));
            this->buffer_z_[i] = static_cast<float>(math::Z(sys.position(i)));
        }
        const std::int32_t block_size(sizeof(float) * sys.size());
        {
            ofs.write(as_bytes(block_size), sizeof(std::int32_t));
            for(std::size_t i=0; i<sys.size(); ++i)
            {
                ofs.write(as_bytes(this->buffer_x_.at(i)), sizeof(float));
            }
            ofs.write(as_bytes(block_size), sizeof(std::int32_t));
        }
        {
            ofs.write(as_bytes(block_size), sizeof(std::int32_t));
            for(std::size_t i=0; i<sys.size(); ++i)
            {
                ofs.write(as_bytes(this->buffer_y_.at(i)), sizeof(float));
            }
            ofs.write(as_bytes(block_size), sizeof(std::int32_t));
        }
        {
            ofs.write(as_bytes(block_size), sizeof(std::int32_t));
            for(std::size_t i=0; i<sys.size(); ++i)
            {
                ofs.write(as_bytes(this->buffer_z_.at(i)), sizeof(float));
            }
            ofs.write(as_bytes(block_size), sizeof(std::int32_t));
        }
    }

    // ------------------------------------------------------------------------
    // write velocity
    {
        std::ofstream ofs(this->vel_name_, std::ios::app | std::ios::binary);
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            this->buffer_x_[i] = static_cast<float>(math::X(sys.velocity(i)));
            this->buffer_y_[i] = static_cast<float>(math::Y(sys.velocity(i)));
            this->buffer_z_[i] = static_cast<float>(math::Z(sys.velocity(i)));
        }
        const std::int32_t block_size(sizeof(float) * sys.size());
        {
            ofs.write(as_bytes(block_size), sizeof(std::int32_t));
            for(std::size_t i=0; i<sys.size(); ++i)
            {
                ofs.write(as_bytes(this->buffer_x_.at(i)), sizeof(float));
            }
            ofs.write(as_bytes(block_size), sizeof(std::int32_t));
        }
        {
            ofs.write(as_bytes(block_size), sizeof(std::int32_t));
            for(std::size_t i=0; i<sys.size(); ++i)
            {
                ofs.write(as_bytes(this->buffer_y_.at(i)), sizeof(float));
            }
            ofs.write(as_bytes(block_size), sizeof(std::int32_t));
        }
        {
            ofs.write(as_bytes(block_size), sizeof(std::int32_t));
            for(std::size_t i=0; i<sys.size(); ++i)
            {
                ofs.write(as_bytes(this->buffer_z_.at(i)), sizeof(float));
            }
            ofs.write(as_bytes(block_size), sizeof(std::int32_t));
        }
    }

    // ------------------------------------------------------------------------
    // write energy
    {
        std::ofstream ofs(this->ene_name_, std::ios::app);
        // if the width exceeds, operator<<(std::ostream, std::string) ignores
        // ostream::width and outputs whole string.
        ofs << std::setw(11) << std::left << std::to_string(step) << ' ';
        ofs << ff.dump_energy(sys) << ' ';
        ofs << std::setw(14) << std::right << this->calc_kinetic_energy(sys) << '\n';
        ofs.close();

        // XXX consider introducing template argument to remove this if-branching
        //     at the compile time
        if(this->output_progress_)
        {
            std::cerr << progress_bar_.format(step);
            if(step == progress_bar_.total()){std::cerr << std::endl;}
        }
    }
    return ;
}

} // mjolnir
#endif // MJOLNIR_CORE_DCD_OBSERVER_HPP
