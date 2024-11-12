/*
 * eta_progress_bar.h
 *
 * A custom ProgressBar class to display a progress bar with time estimation
 *
 * Author: clemens@nevrome.de
 *
 * Copied from https://github.com/kforner/rcpp_progress/blob/master/inst/examples/RcppProgressETA/src/eta_progress_bar.hpp with modifications by NHG
 *
 */
#ifndef _RcppProgress_ETA_PROGRESS_BAR_HPP
#define _RcppProgress_ETA_PROGRESS_BAR_HPP

#include <R_ext/Print.h>
#include <ctime>
#include <stdio.h>
#include <sstream>
#include <string.h>

#include <progress.hpp>
#include "progress_bar.hpp"

// for unices only
#if !defined(WIN32) && !defined(__WIN32) && !defined(__WIN32__)
#include <Rinterface.h>
#endif

class ETAProgressBar: public ProgressBar{
  public: // ====== LIFECYCLE =====

    /**
     * Main constructor
     */
    ETAProgressBar()  {
      _max_ticks = 50;
      _finalized = false;
      _timer_flag = true;
    }

    ~ETAProgressBar() {
    }

    public: // ===== main methods =====

      void display() {
        REprintf("0%%   10   20   30   40   50   60   70   80   90   100%%\n");
        REprintf("[----|----|----|----|----|----|----|----|----|----|\n");
        flush_console();
      }

      // update display
      void update(float progress) {

        // stop if already finalized
        if (_finalized) return;

        // measure current time
        time(&current_time);

        // start time measurement when update() is called the first time
        if (_timer_flag) {
          _timer_flag = false;
          // measure start time
          time_at_start = current_time;

          time_at_last_refresh = current_time;

          progress_at_last_refresh = progress;

          _num_ticks = _compute_nb_ticks(progress);

          time_string = "calculating...";

          // create progress bar string
          std::string progress_bar_string = _current_ticks_display(_num_ticks);

          // merge progress bar and time string
          std::stringstream strs;
          strs << "|" << progress_bar_string << "| ETA: " << time_string;
          std::string temp_str = strs.str();
          char const* char_type = temp_str.c_str();

          // print: remove old and replace with new
          REprintf("\r");
          REprintf("%s", char_type);
        } else {

          double time_since_start = std::difftime(current_time, time_at_start);

          if (progress != 1) {
            // ensure overwriting of old time info
            int empty_length = time_string.length();

            int _num_ticks_current = _compute_nb_ticks(progress);

            bool update_bar = (_num_ticks_current != _num_ticks);

            _num_ticks = _num_ticks_current;

            if (progress > 0 && time_since_start > 1) {
              double time_since_last_refresh = std::difftime(current_time, time_at_last_refresh);

              if (time_since_last_refresh >= .5) {
                update_bar = true;

                double progress_since_last_refresh = progress - progress_at_last_refresh;

                double total_rate = progress / time_since_start;

                if (progress_since_last_refresh == 0) {
                  progress_since_last_refresh = .0000001;
                }

                double current_rate = progress_since_last_refresh / time_since_last_refresh;

                //alpha weights average rate against recent recent (current) rate;
                //alpha = 1 => estimate based on on total_rate (treats as constant)
                //alpha = 0 => estimate based on recent rate (high fluctuation)
                double alpha = .8;

                double eta = (1 - progress) * (alpha / total_rate + (1 - alpha) / current_rate);

                // convert seconds to time string
                time_string = "~";
                time_string += _time_to_string(eta);

                time_at_last_refresh = current_time;
                progress_at_last_refresh = progress;
              }
            }

            if (update_bar) {
              // create progress bar string
              std::string progress_bar_string = _current_ticks_display(_num_ticks);

              std::string empty_space = std::string(std::fdim(empty_length, time_string.length()), ' ');

              // merge progress bar and time string
              std::stringstream strs;
              strs << "|" << progress_bar_string << "| ETA: " << time_string << empty_space;
              std::string temp_str = strs.str();
              char const* char_type = temp_str.c_str();

              // print: remove old and replace with new
              REprintf("\r");
              REprintf("%s", char_type);
            }

          } else {
            // ensure overwriting of old time info
            int empty_length = time_string.length();

            // finalize display when ready

            // convert seconds to time string
            std::string time_string = _time_to_string(time_since_start);

            std::string empty_space = std::string(std::fdim(empty_length, time_string.length()), ' ');

            // create progress bar string
            _num_ticks = _compute_nb_ticks(progress);

            std::string progress_bar_string = _current_ticks_display(_num_ticks);

            // merge progress bar and time string
            std::stringstream strs;
            strs << "|" << progress_bar_string << "| " << "Elapsed: " << time_string << empty_space;

            std::string temp_str = strs.str();
            char const* char_type = temp_str.c_str();

            // print: remove old and replace with new
            REprintf("\r");
            REprintf("%s", char_type);

            _finalize_display();
          }
        }
      }

      void end_display() {
        update(1);
      }

      protected: // ==== other instance methods =====

        // convert double with seconds to time string
        std::string _time_to_string(double seconds) {

          int time = (int) seconds;

          int hour = 0;
          int min = 0;
          int sec = 0;

          hour = time / 3600;
          time = time % 3600;
          min = time / 60;
          time = time % 60;
          sec = time;

          std::stringstream time_strs;
          if (hour != 0) time_strs << hour << "h ";
          if (min != 0) time_strs << min << "min ";
          if (sec != 0 || (hour == 0 && min == 0)) time_strs << sec << "s ";
          std::string time_str = time_strs.str();

          return time_str;
        }

        // update the ticks display corresponding to progress
        std::string _current_ticks_display(int nb) {

          std::stringstream ticks_strs;
          for (int i = 0; i < (_max_ticks - 1); ++i) {
            if (i < nb) {
              ticks_strs << "=";
            } else {
              ticks_strs << " ";
            }
          }
          std::string tick_space_string = ticks_strs.str();

          return tick_space_string;
        }

        // finalize
        void _finalize_display() {
          if (_finalized) return;

          REprintf("\n");
          flush_console();
          _finalized = true;
        }

        // compute number of ticks according to progress
        int _compute_nb_ticks(float progress) {
          return int(progress * _max_ticks);
        }

        // N.B: does nothing on windows
        void flush_console() {
#if !defined(WIN32) && !defined(__WIN32) && !defined(__WIN32__)
          R_FlushConsole();
#endif
        }

        private: // ===== INSTANCE VARIABLES ====
          int _max_ticks;   		// the total number of ticks to print
          int _num_ticks;
          bool _finalized;
          bool _timer_flag;
          time_t time_at_start, current_time, time_at_last_refresh;
          float progress_at_last_refresh;
          std::string time_string;

};

#endif