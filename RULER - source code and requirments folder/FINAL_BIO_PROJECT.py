import cv2
import numpy as np
import pandas as pd
import math
import os
from tkinter import Tk, filedialog
from skimage import io, color, filters, morphology, measure, util
from skimage.morphology import skeletonize
from scipy.ndimage import binary_fill_holes
from scipy.spatial.distance import cdist, pdist, squareform

# Threshold Mode hand trace Calibration (Manual Threshold Mode), the best conversion rate 
# in order to match the known sperm length based on different threshold components being different sizes
THRESHOLD_PIXELS_PER_MICRON = 41.2
THRESHOLD_TRIM_RADIUS = 6 

# Auto mode calibration (measuring using component detection + selection + skeletonization)
AUTO_PIXELS_PER_UM = 3.06
MIN_SPERM_LENGTH_PX = 1000
MIN_FRAGMENT_LENGTH_PX = 10 





# Branch adjustment, pruning, and measuring logic (once the components are found and selected by the user, their lengths are measured here)
def find_skeleton_endpoints(skeleton_image: np.ndarray) -> list[tuple[int, int]]:
    """Identify endpoints in the skeleton (pixels with exactly 1 neighbor)."""
    endpoints = []
    img_height, img_width = skeleton_image.shape
    
    # Find endpoints (pixels with exactly 1 neighbor)
    for current_y in range(img_height):
        for current_x in range(img_width):
            if not skeleton_image[current_y, current_x]: continue
            
            neighbor_count = 0
            for offset_y in [-1, 0, 1]:
                for offset_x in [-1, 0, 1]:
                    if offset_y == 0 and offset_x == 0: continue
                    neighbor_y, neighbor_x = current_y + offset_y, current_x + offset_x
                    
                    if 0 <= neighbor_y < img_height and 0 <= neighbor_x < img_width and skeleton_image[neighbor_y, neighbor_x]:
                        neighbor_count += 1
                        
            if neighbor_count == 1:
                endpoints.append((current_y, current_x))
    return endpoints

def trace_and_remove_branch(skeleton_image: np.ndarray, start_point: tuple[int, int], min_branch_length: int) -> bool:
    """Trace a branch from an endpoint and remove it if it's too short."""
    img_height, img_width = skeleton_image.shape
    start_y, start_x = start_point
    current_branch_pixels = [(start_y, start_x)]
    current_pixel = (start_y, start_x)
    
    # Trace branch until we hit a junction or run out
    while len(current_branch_pixels) < min_branch_length:
        current_y, current_x = current_pixel
        next_pixel_in_branch = None
        
        # Search for the next connected pixel not already in our traced list
        for offset_y in [-1, 0, 1]:
            for offset_x in [-1, 0, 1]:
                if offset_y == 0 and offset_x == 0: continue
                
                neighbor_y, neighbor_x = current_y + offset_y, current_x + offset_x
                
                if 0 <= neighbor_y < img_height and 0 <= neighbor_x < img_width and skeleton_image[neighbor_y, neighbor_x]:
                    if (neighbor_y, neighbor_x) not in current_branch_pixels:
                        next_pixel_in_branch = (neighbor_y, neighbor_x)
                        break
            if next_pixel_in_branch: break
        
        if not next_pixel_in_branch: break
        
        # Check if the next pixel is a junction (more than 2 neighbors)
        next_y, next_x = next_pixel_in_branch
        next_neighbor_count = 0
        for offset_y in [-1, 0, 1]:
            for offset_x in [-1, 0, 1]:
                if offset_y == 0 and offset_x == 0: continue
                
                neighbor_y, neighbor_x = next_y + offset_y, next_x + offset_x
                if 0 <= neighbor_y < img_height and 0 <= neighbor_x < img_width and skeleton_image[neighbor_y, neighbor_x]:
                    next_neighbor_count += 1
        
        current_branch_pixels.append(next_pixel_in_branch)
        current_pixel = next_pixel_in_branch
        
        # If we hit a junction (3+ neighbors), stop tracing this branch
        if next_neighbor_count > 2: break
    
    # If the branch traced is shorter than the minimum length, delete it
    if len(current_branch_pixels) < min_branch_length:
        for pixel_y, pixel_x in current_branch_pixels:
            skeleton_image[pixel_y, pixel_x] = False
        return True # Indicate that a change was made
        
    return False

def remove_small_sperm_branches(skeleton_image: np.ndarray, min_branch_length: int = MIN_FRAGMENT_LENGTH_PX) -> np.ndarray:
    """Remove short branches from skeleton, keeping only main backbone."""
    pruned_skeleton = skeleton_image.copy()
    has_changes = True
    
    while has_changes:
        has_changes = False
        branch_endpoints = find_skeleton_endpoints(pruned_skeleton)
        
        # For each endpoint, trace back and remove if branch is short
        for start_point in branch_endpoints:
             if trace_and_remove_branch(pruned_skeleton, start_point, min_branch_length):
                 has_changes = True
                
    return pruned_skeleton


def measure_len_of_sperm_skeleton(skeleton_image: np.ndarray) -> float:
    """Compute length of a 1-pixel-wide skeleton (New Approach)."""
    # Ensure the skeleton is boolean for logical operations
    if skeleton_image.dtype != bool:
        skeleton_image = skeleton_image.astype(bool)

    img_height, img_width = skeleton_image.shape
    total_skeleton_length = 0.0

    # Define 4-connected neighbors (half-set) to count edges exactly once.
    # (Right, Down, Down-Right, Down-Left)
    # Distance is 1.0 for direct neighbors, sqrt(2) for diagonal neighbors.
    half_neighbor_offsets = [
        (0, 1, 1.0),            # Right
        (1, 0, 1.0),            # Down
        (1, 1, math.sqrt(2)),   # Down-Right
        (1, -1, math.sqrt(2)),  # Down-Left
    ]

    # Get coordinates of all skeleton pixels
    skeleton_y_coords, skeleton_x_coords = np.nonzero(skeleton_image)

    # Iterate through every pixel in the skeleton
    for current_y, current_x in zip(skeleton_y_coords, skeleton_x_coords):
        # Check specific neighbors to accumulate length
        for offset_y, offset_x, neighbor_distance in half_neighbor_offsets:
            neighbor_y = current_y + offset_y
            neighbor_x = current_x + offset_x
            
            # Boundary check: Ensure neighbor is within image bounds
            if 0 <= neighbor_y < img_height and 0 <= neighbor_x < img_width:
                # If the neighbor is also part of the skeleton, add the distance
                if skeleton_image[neighbor_y, neighbor_x]:
                    total_skeleton_length += neighbor_distance

    return total_skeleton_length



# --- MAIN CLASS CONTAINING UI AND THRESHOLD MODE AND GUI LOGIC ---

class SpermSizer:
    def __init__(self, image_path):
        self.image_path = image_path
        
        # Load initial image for display mode using OpenCV
        self.img_gray_original = cv2.imread(image_path, cv2.IMREAD_GRAYSCALE)
        if self.img_gray_original is None:
            raise ValueError(f"Image not found or path invalid: {image_path}")
        
        self.inverted = False
        self.img_gray = self.img_gray_original.copy()
        
        self.sperm_cells = []
        self.overlay = cv2.cvtColor(self.img_gray.copy(), cv2.COLOR_GRAY2BGR)
        self.running = True
        
        # Auto measure mode 
        self.auto_overlay = None
        self.auto_length_px = 0.0
        self.auto_length_um = 0.0
        self.auto_message = "Press 'm' to Measure"
        self.auto_segments = [] 
        self.auto_selected_indices = set() 
        self.auto_threshold_adjust = 1.0 
        
        # Auto measure mode Manual Editing Tools
        self.auto_add_line_mode = False
        self.auto_add_pixel_mode = False
        self.auto_remove_mode = False
        self.auto_line_start = None
        self.auto_manual_lines = []  # [(x1, y1, x2, y2, length_px), ...]
        self.auto_manual_pixels = [] # [(x, y), ...]

        # Manual thresholding hand trace approach
        self.threshold_value = 128
        self.threshold_mode = False # Start in Auto Mode by default
        
        # Threshold Mode Tools
        self.pixel_measuring_mode = False
        self.line_mode = False
        self.line_start_point = None
        self.line_preview_end = None
        
        # Storage for Threshold Mode
        self.highlighted_components = []
        self.highlighted_pixels = [] 
        self.highlighted_lines = []   
        self.measured_components = {} 
        self.measured_pixels = {}     
        self.measured_lines = {}      
        
        # Visualization settings
        self.color_index = 0
        self.distinct_colors = [
            [0, 255, 255], [255, 0, 255], [255, 255, 0], [0, 255, 0],
            [0, 0, 255], [255, 0, 0], [0, 165, 255], [203, 192, 255],
            [130, 0, 75], [0, 255, 127], [255, 255, 255], [128, 128, 128],
        ]
        self.selecting = False
        self.selection_start = None
        self.selection_end = None
        self.show_highlights = True
        self.current_threshold_binary = None

    #Auto measure 
    
    def run_auto_measure(self):
        """Runs the Sato/Ridge detection + Skeletonization pipeline."""
        print(f"Running Auto-Measurement pipeline ")#(Sensitivity: {self.auto_threshold_adjust:.2f})...") #NEVIN CHANGE, commented
        self.auto_message = "Calculating..."
        
        try:
            # Prepare image
            img = io.imread(self.image_path)
            img = util.img_as_float(img)
            
            if img.ndim == 3:
                img_gray_auto = color.rgb2gray(img)
            else:
                img_gray_auto = img
            
            if self.inverted:
                img_gray_auto = 1.0 - img_gray_auto
            
            #  Sato Filter
            vessel = filters.sato(img_gray_auto, sigmas=(1, 2, 3), black_ridges=False)
            v_min, v_max = vessel.min(), vessel.max()
            vessel = (vessel - v_min) / (v_max - v_min + 1e-8)
            
            # Threshold (Otsu + Adjustment)
            otsu = filters.threshold_otsu(vessel)
            threshold_val = otsu * self.auto_threshold_adjust
            binary = vessel > threshold_val
            
            # Clean
            # Use same settings as reference to catch EVERYTHING
            binary = morphology.remove_small_objects(binary, min_size=MIN_FRAGMENT_LENGTH_PX)
            binary = morphology.remove_small_holes(binary, area_threshold=50)
            binary = morphology.binary_opening(binary, morphology.disk(1))
            
            # Reset detection state
            self.auto_overlay = cv2.cvtColor(self.img_gray.copy(), cv2.COLOR_GRAY2BGR)
            self.auto_segments = []
            self.auto_selected_indices = set()

            if not binary.any():
                self.update_auto_results()
                self.auto_message = "No structure detected."
                return

            # GLOBAL SKELETONIZATION FOR CONNECTIVITY 
            skel_all = morphology.skeletonize(binary)
            labels = measure.label(skel_all)
            
            # OPTIMIZATION: FILTER SMALL SKELETONS INSTANTLY ---
            component_sizes = np.bincount(labels.ravel())
            valid_labels = np.where(component_sizes > MIN_FRAGMENT_LENGTH_PX)[0]
            valid_labels = valid_labels[valid_labels > 0] # Remove background
            
            max_len = 0
            best_idx = -1
            
            for lab in valid_labels:
                comp_skel = (labels == lab)
                
                # Prune small noise spurs
                comp_pruned = remove_small_sperm_branches(comp_skel, min_branch_length=MIN_FRAGMENT_LENGTH_PX)
                
                # Fallback if pruning kills small segments
                if not comp_pruned.any():
                    comp_pruned = comp_skel
                
                L = measure_len_of_sperm_skeleton(comp_pruned)
                path_skel = comp_pruned
                
                if L >= MIN_FRAGMENT_LENGTH_PX:
                    seg_id = len(self.auto_segments)
                    ys, xs = np.nonzero(path_skel)
                    
                    # Store data
                    self.auto_segments.append({
                        'id': seg_id,
                        'pixels': (ys, xs),
                        'length': L,
                        'path_mask': path_skel,
                        'display_mask': morphology.dilation(path_skel, morphology.disk(1))
                    })
                    
                    if L > max_len:
                        max_len = L
                        best_idx = seg_id
            
            # Select the longest segment by default
            if best_idx != -1:
                self.auto_selected_indices.add(best_idx)
                
            self.update_auto_results()
            self.draw_auto_overlay()

        except Exception as e:
            print(f"Auto-measure error: {e}")
            self.auto_message = "Error in auto-measure"

    def update_auto_results(self):
        #  Sum Auto Segments
        auto_px = sum(self.auto_segments[i]['length'] for i in self.auto_selected_indices)
        
        #  Sum Manual Lines
        manual_line_px = sum(line[4] for line in self.auto_manual_lines)
        
        # Sum Manual Pixels (Approximate 1 pixel length per dot)
        manual_dot_px = len(self.auto_manual_pixels)
        
        total_px = auto_px + manual_line_px + manual_dot_px
        self.auto_length_px = total_px
        self.auto_length_um = total_px / AUTO_PIXELS_PER_UM
        
        count_auto = len(self.auto_selected_indices)
        count_lines = len(self.auto_manual_lines)
        count_dots = len(self.auto_manual_pixels)
        
        if total_px == 0:
             self.auto_message = "No measurement. Click 'm' or add lines."
        else:
             self.auto_message = f"L: {self.auto_length_um:.2f} um (Segs:{count_auto}, Lines:{count_lines}, Dots:{count_dots})"

    def draw_auto_overlay(self):
        # Start with fresh base image
        base = cv2.cvtColor(self.img_gray.copy(), cv2.COLOR_GRAY2BGR)
        
        #  Draw Auto Segments
        for seg in self.auto_segments:
            idx = seg['id']
            if idx in self.auto_selected_indices:
                color = (0, 255, 0) # Green (Selected)
            else:
                color = (255, 0, 0) # Blue (Detected/Ignored)
            
            # Use the display mask (dilated skeleton) to draw FILLED contours
            mask_uint8 = (seg['display_mask'].astype(np.uint8)) * 255
            # Use RETR_TREE to capture holes (loops)
            contours, _ = cv2.findContours(mask_uint8, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
            
            # Draw outlines (thickness 1) to reveal loops/holes inside
            cv2.drawContours(base, contours, -1, color, 1)
            
        # Draw Manual Lines (Yellow)
        for x1, y1, x2, y2, _ in self.auto_manual_lines:
            cv2.line(base, (x1, y1), (x2, y2), (0, 255, 255), 2)
            
        # Draw Manual Pixels (Cyan)
        for x, y in self.auto_manual_pixels:
            cv2.circle(base, (x, y), 2, (255, 255, 0), -1)
            
        # Draw Line Preview
        if self.auto_add_line_mode and self.auto_line_start:
            pass 

        self.auto_overlay = base

    def check_click_auto_mode(self, x, y):
        """Handle clicks in Auto Mode based on current tool state."""
        
        # REMOVE MODE
        if self.auto_remove_mode:
            # Check Manual Pixels
            for i, (px, py) in enumerate(self.auto_manual_pixels):
                if abs(px - x) < 5 and abs(py - y) < 5:
                    self.auto_manual_pixels.pop(i)
                    self.update_auto_results()
                    self.draw_auto_overlay()
                    return

            # Check Manual Lines (Distance to segment)
            for i, (x1, y1, x2, y2, _) in enumerate(self.auto_manual_lines):
                px = x2 - x1; py = y2 - y1
                norm = px*px + py*py
                u =  ((x - x1) * px + (y - y1) * py) / float(norm) if norm != 0 else 0
                if u > 1: u = 1
                elif u < 0: u = 0
                dx = x1 + u * px; dy = y1 + u * py
                dist = math.sqrt((x - dx)**2 + (y - dy)**2)
                
                if dist < 5:
                    self.auto_manual_lines.pop(i)
                    self.update_auto_results()
                    self.draw_auto_overlay()
                    return

            #  Check Auto Segments (Deselect)
            clicked_id = -1
            
            # Check inside mask first
            for seg in self.auto_segments:
                mask = seg['display_mask']
                if 0 <= y < mask.shape[0] and 0 <= x < mask.shape[1] and mask[y, x]:
                    clicked_id = seg['id']
                    break
            
            # Check distance fallback
            if clicked_id == -1:
                min_dist = 15
                for seg in self.auto_segments:
                    ys, xs = seg['pixels']
                    dists = np.sqrt((xs - x)**2 + (ys - y)**2)
                    cur_min = dists.min()
                    if cur_min < min_dist:
                        min_dist = cur_min
                        clicked_id = seg['id']
            
            if clicked_id != -1 and clicked_id in self.auto_selected_indices:
                self.auto_selected_indices.remove(clicked_id)
                self.update_auto_results()
                self.draw_auto_overlay()
            return

        #  ADD LINE MODE 
        if self.auto_add_line_mode:
            if self.auto_line_start is None:
                self.auto_line_start = (x, y)
                print(f"Line start: {x},{y}")
            else:
                x1, y1 = self.auto_line_start
                dist_px = math.sqrt((x - x1)**2 + (y - y1)**2)
                self.auto_manual_lines.append((x1, y1, x, y, dist_px))
                self.auto_line_start = None
                self.update_auto_results()
                self.draw_auto_overlay()
            return

        # ADD PIXEL MODE 
        if self.auto_add_pixel_mode:
            self.auto_manual_pixels.append((x, y))
            self.update_auto_results()
            self.draw_auto_overlay()
            return

        # DEFAULT: SELECT/DESELECT AUTO SEGMENTS
        clicked_id = -1
        # Check inside mask
        for seg in self.auto_segments:
            mask = seg['display_mask']
            if 0 <= y < mask.shape[0] and 0 <= x < mask.shape[1] and mask[y, x]:
                clicked_id = seg['id']
                break
        
        # Fallback distance check to skeleton (for thin parts)
        if clicked_id == -1:
            min_dist = 15
            for seg in self.auto_segments:
                ys, xs = seg['pixels']
                dists = np.sqrt((xs - x)**2 + (ys - y)**2)
                cur_min = dists.min()
                if cur_min < min_dist:
                    min_dist = cur_min
                    clicked_id = seg['id']
        
        if clicked_id != -1:
            if clicked_id in self.auto_selected_indices:
                self.auto_selected_indices.remove(clicked_id)
            else:
                self.auto_selected_indices.add(clicked_id)
            self.update_auto_results()
            self.draw_auto_overlay()

    #  Manual  EXISTING THRESHOLD MODE LOGIC

    def ask_invert(self):
        pass 
    
    def toggle_inversion(self):
        self.inverted = not self.inverted
        if self.inverted:
            self.img_gray = cv2.bitwise_not(self.img_gray_original)
        else:
            self.img_gray = self.img_gray_original.copy()
        
        if not self.threshold_mode:
            self.auto_overlay = None
            self.auto_segments = []
            self.auto_selected_indices = set()
            self.auto_manual_lines = []
            self.auto_manual_pixels = []
            self.auto_message = "Image inverted. Press 'm' to Measure."
        
        self.update_overlay()
        self.highlighted_components = []
        self.highlighted_pixels = []
        self.measured_pixels = {}
    
    def update_overlay(self):
        self.overlay = cv2.cvtColor(self.img_gray.copy(), cv2.COLOR_GRAY2BGR)
        if self.show_highlights:
            self.apply_highlights_to_image(self.overlay)
        
    def apply_highlights_to_image(self, target_image):
        """Applies manual highlights (components, pixels, lines) to ANY given image."""
        for mask, thresh_val, color_idx in self.highlighted_components:
            color = self.distinct_colors[color_idx % len(self.distinct_colors)]
            highlight_overlay = target_image.copy()
            highlight_overlay[mask > 0] = color
            cv2.addWeighted(highlight_overlay, 0.6, target_image, 0.4, 0, target_image)
        
        for x, y, color_idx in self.highlighted_pixels:
            color = self.distinct_colors[color_idx % len(self.distinct_colors)]
            cv2.circle(target_image, (x, y), 3, color, -1)
        
        for x1, y1, x2, y2, color_idx in self.highlighted_lines:
            color = self.distinct_colors[color_idx % len(self.distinct_colors)]
            cv2.line(target_image, (x1, y1), (x2, y2), color, 2)
    
    def load_new_image(self):
        print("Select a new image file...")
        new_img_path = filedialog.askopenfilename(
            title="Choose sperm image",
            filetypes=[("Images", "*.jpg *.jpeg *.png *.bmp *.tif *.tiff"), ("All files", "*")]
        )
        
        if new_img_path:
            if len(self.sperm_cells) > 0:
                self.export_results()
            
            self.image_path = new_img_path
            self.img_gray_original = cv2.imread(new_img_path, cv2.IMREAD_GRAYSCALE)
            
            if self.img_gray_original is None:
                print(f"Error: Could not load image from {new_img_path}")
                return
            
            if self.inverted:
                self.img_gray = cv2.bitwise_not(self.img_gray_original)
            else:
                self.img_gray = self.img_gray_original.copy()
            
            self.sperm_cells = []
            self.highlighted_components = []
            self.highlighted_pixels = []
            self.highlighted_lines = []
            self.measured_pixels = {}
            self.measured_lines = {}
            self.measured_components = {}
            self.threshold_mode = False # Reset to Auto Mode on load
            self.pixel_measuring_mode = False
            self.line_mode = False
            self.line_start_point = None
            self.line_preview_end = None
            
            # Reset Auto Measure State
            self.auto_overlay = None
            self.auto_segments = []
            self.auto_selected_indices = set()
            self.auto_manual_lines = []
            self.auto_manual_pixels = []
            self.auto_message = "Press 'm' to Measure"
            
            self.update_overlay()
            print(f"Loaded new image: {new_img_path}")

    def create_instructions_image(self):
        img_height = 1000
        img_width = 1200
        instructions_img = np.zeros((img_height, img_width, 3), dtype=np.uint8)
        instructions_img[:] = (40, 40, 40)
        
        instructions = [
            "SPERM SIZER INSTRUCTIONS",
            "",
            "GENERAL CONTROLS:",
            "  'i' - Toggle image inversion",
            "  'l' - Load new image",
            "  't' - Toggle Mode (Auto vs Manual Threshold)",
            "  'h' - Toggle highlight visibility",
            "  'p' - Print color summary to console",
            "  'e' - Export color summary to CSV",
            "  ESC - Exit application",
            "",
        ]
        
        if self.inverted:
            instructions.append("STATUS: Image currently INVERTED")
        
        if self.threshold_mode:
            pixel_mode_text = "PIXEL MEASURING MODE: ON" if self.pixel_measuring_mode else "Pixel Measuring Mode: OFF"
            line_mode_text = "LINE MODE: ON" if self.line_mode else "Line Mode: OFF"
            instructions.extend([
                "--- MANUAL THRESHOLD MODE (Active) ---",
                f"  Current threshold: {self.threshold_value}",
                "  '+'/'-' or UP/DOWN - Adjust threshold",
                "  's' - Save thresholded image",
                f"  Current color: #{self.color_index + 1}",
                "  '1'-'9','0' - Select color",
                f"  {pixel_mode_text}",
                "  'u' - Toggle pixel measuring mode",
                f"  {line_mode_text}",
                "  'k' - Toggle line drawing mode",
                "  In line mode: Click start, then end point to draw line",
                "  Click sperm - Highlight single component/pixel",
                "  Click & Drag - Highlight all components in region",
                "  Click highlighted - Remove highlight",
                "  'c' - Clear all highlights",
            ])
        else:
            add_line_txt = "ON" if self.auto_add_line_mode else "OFF"
            add_pix_txt = "ON" if self.auto_add_pixel_mode else "OFF"
            remove_txt = "ON" if self.auto_remove_mode else "OFF"
            
            instructions.extend([
                "--- AUTO MEASURE MODE (Active) ---",
                "  (Press 't' to switch to Manual Threshold Mode)",
                "",
                f"  STATUS: {self.auto_message}",
                f"  SENSITIVITY: {self.auto_threshold_adjust:.2f} (Default 1.00)",
                "",
                "  CONTROLS:",
                "  'm' - Run Auto-Measurement Pipeline",
                "  '[' / ']' - Decrease/Increase Sensitivity",
                "  CLICK on segments to Select/Deselect",
                "",
                "  EDITING TOOLS:",
                f"  'n' - Add Line Mode ({add_line_txt}) - Click 2 points",
                f"  'b' - Add Pixel Mode ({add_pix_txt}) - Click point",
                f"  'r' - Remove Mode ({remove_txt}) - Click item to remove",
                "",
                "  LEGEND:",
                "   GREEN: Selected Auto-Segment",
                "   BLUE: Ignored Auto-Segment",
                "   YELLOW: Manual Line",
                "   CYAN: Manual Pixel",
            ])
            
        y_pos = 30
        for line in instructions:
            color = (255, 255, 255)
            if "Active" in line: color = (0, 255, 0)
            elif "===" in line: color = (0, 255, 255)
            
            cv2.putText(instructions_img, line, (20, y_pos), 
                       cv2.FONT_HERSHEY_SIMPLEX, 0.6, color, 1, cv2.LINE_AA)
            y_pos += 25
        
        # Add measurements section on the right side
        measurements_x = 500
        measurements_y = 30
        
        # Display logic based on mode
        if not self.threshold_mode:
             cv2.putText(instructions_img, "    AUTO RESULTS    ", (measurements_x, measurements_y),
                        cv2.FONT_HERSHEY_SIMPLEX, 0.7, (0, 255, 255), 1, cv2.LINE_AA)
             measurements_y += 40
             if self.auto_length_px > 0:
                 cv2.putText(instructions_img, f"Length: {self.auto_length_um:.2f} um", (measurements_x, measurements_y),
                            cv2.FONT_HERSHEY_SIMPLEX, 0.8, (0, 255, 0), 2, cv2.LINE_AA)
                 measurements_y += 30
                 cv2.putText(instructions_img, f"Pixels: {self.auto_length_px:.1f} px", (measurements_x, measurements_y),
                            cv2.FONT_HERSHEY_SIMPLEX, 0.6, (200, 200, 200), 1, cv2.LINE_AA)
                 measurements_y += 30
                 cv2.putText(instructions_img, f"Segs: {len(self.auto_selected_indices)} | Lines: {len(self.auto_manual_lines)} | Dots: {len(self.auto_manual_pixels)}", 
                             (measurements_x, measurements_y), cv2.FONT_HERSHEY_SIMPLEX, 0.6, (200, 200, 200), 1, cv2.LINE_AA)
             else:
                 cv2.putText(instructions_img, "No data yet. Press 'm'.", (measurements_x, measurements_y),
                            cv2.FONT_HERSHEY_SIMPLEX, 0.6, (200, 200, 200), 1, cv2.LINE_AA)
        
        # In both modes, show manual measurements if they exist
        if len(self.measured_components) > 0 or len(self.measured_pixels) > 0:
            # Adjust Y pos if we already printed auto results
            if not self.threshold_mode: measurements_y += 40
            
            cv2.putText(instructions_img, "   MANUAL MODE DATA   ", (measurements_x, measurements_y),
                        cv2.FONT_HERSHEY_SIMPLEX, 0.7, (0, 255, 255), 1, cv2.LINE_AA)
            measurements_y += 40
            
            # Add color sums section
            color_sums = self.get_color_sums()
            if color_sums:
                 for color_idx in sorted(color_sums.keys()):
                    data = color_sums[color_idx]
                    color_rgb = self.distinct_colors[color_idx % len(self.distinct_colors)]
                    # Create a small color swatch
                    swatch_y = measurements_y - 5
                    cv2.rectangle(instructions_img, (measurements_x, swatch_y), 
                                (measurements_x + 20, swatch_y + 15), tuple(map(int, color_rgb)), -1)
                    summary_text = f"Color #{color_idx + 1}: {data['count']} cells, {data['sum_um']:.2f} um"
                    cv2.putText(instructions_img, summary_text, (measurements_x + 25, measurements_y + 10),
                               cv2.FONT_HERSHEY_SIMPLEX, 0.5, (200, 200, 200), 1, cv2.LINE_AA)
                    measurements_y += 25
        
        # DRAW COLOR PALETTE CHART To determine the color of the component 
        if self.threshold_mode:
            palette_y = img_height - 100
            start_x = 20
            box_size = 40
            spacing = 10
            
            cv2.putText(instructions_img, "COLOR PALETTE (Keys 1-9, 0):", (start_x, palette_y - 15), 
                        cv2.FONT_HERSHEY_SIMPLEX, 0.6, (255, 255, 255), 1, cv2.LINE_AA)

            for i in range(10):
                color = self.distinct_colors[i]
                x = start_x + i * (box_size + spacing)
                cv2.rectangle(instructions_img, (x, palette_y), (x + box_size, palette_y + box_size), 
                            tuple(map(int, color)), -1)
                key_num = (i + 1) % 10
                text = str(key_num)
                cv2.putText(instructions_img, text, (x + 12, palette_y + 28), 
                            cv2.FONT_HERSHEY_SIMPLEX, 0.7, (0, 0, 0), 3) 
                cv2.putText(instructions_img, text, (x + 12, palette_y + 28), 
                            cv2.FONT_HERSHEY_SIMPLEX, 0.7, (255, 255, 255), 1) 
                
                if i == self.color_index:
                    cv2.rectangle(instructions_img, (x - 3, palette_y - 3), 
                                (x + box_size + 3, palette_y + box_size + 3), 
                                (255, 255, 255), 3)

        return instructions_img

    def click_event(self, event, x, y, flags, param):
        if self.threshold_mode:
            if event == cv2.EVENT_LBUTTONDOWN:
                self.selecting = True
                self.selection_start = (x, y)
                self.selection_end = (x, y)
            
            elif event == cv2.EVENT_MOUSEMOVE:
                if self.selecting:
                    self.selection_end = (x, y)
                elif self.line_mode and self.line_start_point is not None:
                    self.line_preview_end = (x, y)
            
            elif event == cv2.EVENT_LBUTTONUP:
                self.selecting = False
                self.selection_end = (x, y)
                start_x, start_y = self.selection_start
                
                if abs(x - start_x) < 5 and abs(y - start_y) < 5:
                    if self.check_click_on_line(x, y):
                        return

                    if self.line_mode:
                        if self.line_start_point is None:
                            self.line_start_point = (x, y)
                            self.line_preview_end = (x, y)
                            print(f"Line start set: ({x}, {y})")
                        else:
                            self.highlight_line(self.line_start_point[0], self.line_start_point[1], x, y)
                            self.line_start_point = None
                            self.line_preview_end = None
                    elif self.pixel_measuring_mode:
                        pixel_clicked = False
                        for i, (px, py, color_idx) in enumerate(self.highlighted_pixels):
                            if (px, py) == (x, y):
                                self.highlighted_pixels.pop(i)
                                if (x,y) in self.measured_pixels: del self.measured_pixels[(x,y)]
                                pixel_clicked = True
                                self.update_overlay()
                                break
                        if not pixel_clicked:
                            self.highlight_pixel(x, y)
                    else:
                        self.highlight_connected_component(x, y)
                else:
                    if not self.pixel_measuring_mode and not self.line_mode:
                        self.highlight_components_in_region(self.selection_start, self.selection_end)
                
                self.selection_start = None
                self.selection_end = None
        else:
            # Auto Mode Click Logic
            # Update preview for line mode
            if self.auto_add_line_mode and self.auto_line_start and event == cv2.EVENT_MOUSEMOVE:
                pass 

            if event == cv2.EVENT_LBUTTONUP:
                self.check_click_auto_mode(x, y)

    def generate_threshold_preview(self):
        gray = self.img_gray.copy()
        gray = cv2.GaussianBlur(gray, (5,5), 0)
        clahe = cv2.createCLAHE(clipLimit=3.0, tileGridSize=(8,8))
        gray = clahe.apply(gray)
        
        _, binary = cv2.threshold(gray, self.threshold_value, 255, cv2.THRESH_BINARY_INV)
        
        kernel = np.ones((3,3), np.uint8)
        binary = cv2.morphologyEx(binary, cv2.MORPH_CLOSE, kernel, iterations=2)
        binary = cv2.morphologyEx(binary, cv2.MORPH_OPEN, kernel, iterations=1)
        
        self.current_threshold_binary = binary
        preview = cv2.cvtColor(binary, cv2.COLOR_GRAY2BGR)
        
        if len(self.highlighted_components) > 0 or len(self.highlighted_pixels) > 0 or len(self.highlighted_lines) > 0 or (self.line_mode and self.line_start_point is not None):
            gray_preview = cv2.cvtColor(cv2.cvtColor(preview, cv2.COLOR_BGR2GRAY), cv2.COLOR_GRAY2BGR)
            colored_highlight = gray_preview.copy()
            
            for mask, thresh_val, color_idx in self.highlighted_components:
                color = self.distinct_colors[color_idx % len(self.distinct_colors)]
                colored_highlight[mask > 0] = color
            for x, y, color_idx in self.highlighted_pixels:
                color = self.distinct_colors[color_idx % len(self.distinct_colors)]
                cv2.circle(colored_highlight, (x, y), 3, color, -1)
            for x1, y1, x2, y2, color_idx in self.highlighted_lines:
                color = self.distinct_colors[color_idx % len(self.distinct_colors)]
                cv2.line(colored_highlight, (x1, y1), (x2, y2), color, 2)
            
            if self.line_mode and self.line_start_point is not None and self.line_preview_end is not None:
                color = self.distinct_colors[self.color_index % len(self.distinct_colors)]
                cv2.line(colored_highlight, self.line_start_point, self.line_preview_end, color, 2)
            
            preview = colored_highlight
        
        cv2.putText(preview, f"Threshold: {self.threshold_value}", 
                    (10, preview.shape[0] - 20), 
                    cv2.FONT_HERSHEY_SIMPLEX, 1, (0, 255, 0), 2)
        return preview
    
    def highlight_components_in_region(self, start, end):
        x1, y1 = start
        x2, y2 = end
        min_x, max_x = min(x1, x2), max(x1, x2)
        min_y, max_y = min(y1, y2), max(y1, y2)
        
        if self.current_threshold_binary is None: return

        num_labels, labels = cv2.connectedComponents(self.current_threshold_binary)
        roi_labels = labels[min_y:max_y+1, min_x:max_x+1]
        unique_labels = np.unique(roi_labels)
        
        for label in unique_labels:
            if label == 0: continue
            
            component_mask = (labels == label).astype(np.uint8) * 255
            if np.any(roi_labels == label):
                already = False
                for i, (mask, _, _) in enumerate(self.highlighted_components):
                    if np.array_equal(mask, component_mask):
                        already = True
                        break
                
                if not already:
                    self.highlighted_components.append((component_mask, self.threshold_value, self.color_index))
                    comp_id = hash(component_mask.tobytes())
                    _, l_um = self.measure_component_length(component_mask)
                    self.measured_components[comp_id] = l_um

    def measure_component_length(self, component_mask):
        binary = component_mask.astype(np.uint8)
        kernel = np.ones((3,3), np.uint8)
        binary = cv2.morphologyEx(binary, cv2.MORPH_CLOSE, kernel, iterations=2)
        binary = cv2.morphologyEx(binary, cv2.MORPH_OPEN, kernel, iterations=1)
        filled = binary_fill_holes(binary > 0)
        binary = (filled * 255).astype(np.uint8)
        binary = self.smooth_line(binary)
        skeleton = self.skeletonize_cell(binary)
        skeleton = self.remove_skeleton_branches(skeleton)
        
        ys, xs = np.nonzero(skeleton)
        if len(xs) < 2: return 0, 0
        
        coords = np.stack([xs, ys], axis=1)
        if len(coords) > 1:
            dist_matrix = squareform(pdist(coords))
            max_dist_idx = np.unravel_index(dist_matrix.argmax(), dist_matrix.shape)
            endpoints = [coords[max_dist_idx[0]], coords[max_dist_idx[1]]]
            skeleton = self.trim_endpoints(skeleton, endpoints, THRESHOLD_TRIM_RADIUS)
        
        return self.measure_length(skeleton, THRESHOLD_PIXELS_PER_MICRON)
    
    def remove_skeleton_branches(self, skeleton, min_branch_length=5):
        skel = skeleton.copy()
        kernel = np.array([[1, 1, 1], [1, 0, 1], [1, 1, 1]], dtype=np.uint8)
        neighbors = cv2.filter2D(skel.astype(np.uint8), -1, kernel) * skel
        
        branch_points = (neighbors > 2) & (skel > 0)
        endpoints = (neighbors == 1) & (skel > 0)
        endpoint_coords = np.argwhere(endpoints)
        
        for ep in endpoint_coords:
            visited = set()
            queue = [tuple(ep)]
            path_length = 0
            hit_branch = False
            while queue and path_length < min_branch_length:
                current = queue.pop(0)
                if current in visited: continue
                visited.add(current)
                path_length += 1
                y, x = current
                if branch_points[y, x]:
                    hit_branch = True
                    break
                for dy in [-1, 0, 1]:
                    for dx in [-1, 0, 1]:
                        if dy == 0 and dx == 0: continue
                        ny, nx = y + dy, x + dx
                        if (0 <= ny < skel.shape[0] and 0 <= nx < skel.shape[1] and
                            skel[ny, nx] > 0 and (ny, nx) not in visited):
                            queue.append((ny, nx))
            if not hit_branch and path_length < min_branch_length:
                for y, x in visited:
                    skel[y, x] = 0
        return skel

    def get_color_sums(self):
        color_sums = {}
        for mask, thresh_val, color_idx in self.highlighted_components:
            component_id = hash(mask.tobytes())
            if component_id in self.measured_components:
                length_um = self.measured_components[component_id]
                if color_idx not in color_sums:
                    color_sums[color_idx] = {'sum_um': 0.0, 'count': 0, 'sum_px': 0.0, 'pixel_count': 0}
                color_sums[color_idx]['sum_um'] += length_um
                color_sums[color_idx]['count'] += 1
                color_sums[color_idx]['sum_px'] += length_um * THRESHOLD_PIXELS_PER_MICRON
        
        for x, y, color_idx in self.highlighted_pixels:
            key = (x, y)
            if key in self.measured_pixels:
                length_um = self.measured_pixels[key]
                if color_idx not in color_sums:
                    color_sums[color_idx] = {'sum_um': 0.0, 'count': 0, 'sum_px': 0.0, 'pixel_count': 0}
                color_sums[color_idx]['sum_um'] += length_um
                color_sums[color_idx]['pixel_count'] += 1
                color_sums[color_idx]['sum_px'] += 1.0
        
        for x1, y1, x2, y2, color_idx in self.highlighted_lines:
            key = (x1, y1, x2, y2)
            if key in self.measured_lines:
                length_um = self.measured_lines[key]
                if color_idx not in color_sums:
                    color_sums[color_idx] = {'sum_um': 0.0, 'count': 0, 'sum_px': 0.0, 'pixel_count': 0}
                color_sums[color_idx]['sum_um'] += length_um
                color_sums[color_idx]['count'] += 1 
                color_sums[color_idx]['sum_px'] += length_um * THRESHOLD_PIXELS_PER_MICRON
        return color_sums
    
    def print_color_sums(self):
        color_sums = self.get_color_sums()
        if not color_sums:
            print("No measurements found.")
            return
        print("\n COLOR SUMMARY ")
        for color_idx in sorted(color_sums.keys()):
            data = color_sums[color_idx]
            print(f"Color #{color_idx + 1}: {data['sum_um']:.2f} um")
        print("                    \n")
    
    def export_color_sums(self):
        color_sums = self.get_color_sums()
        if not color_sums: return
        import os
        base = os.path.splitext(os.path.basename(self.image_path))[0]
        out = filedialog.asksaveasfilename(initialfile=f"{base}_colors.csv", defaultextension=".csv")
        if out:
            rows = []
            for idx in sorted(color_sums.keys()):
                d = color_sums[idx]
                rows.append({
                    'Color_Index': idx+1, 
                    'Count': d['count'], 
                    'Pixels': d['pixel_count'], 
                    'Length_um': d['sum_um']
                })
            pd.DataFrame(rows).to_csv(out, index=False)
            print(f"Saved to {out}")

    def check_click_on_highlighted_component(self, x, y):
        for i, (mask, _, _) in enumerate(self.highlighted_components):
            if mask[y, x] > 0:
                cid = hash(mask.tobytes())
                self.highlighted_components.pop(i)
                if cid in self.measured_components: del self.measured_components[cid]
                self.update_overlay()
                return True
        return False
    
    def check_click_on_line(self, x, y, threshold=5):
        for i, (x1, y1, x2, y2, _) in enumerate(self.highlighted_lines):
            A = x - x1; B = y - y1; C = x2 - x1; D = y2 - y1
            dot = A*C + B*D
            len_sq = C*C + D*D
            param = -1
            if len_sq != 0: param = dot / len_sq
            
            if param < 0: xx, yy = x1, y1
            elif param > 1: xx, yy = x2, y2
            else: xx, yy = x1 + param*C, y1 + param*D
            
            dist = np.sqrt((x-xx)**2 + (y-yy)**2)
            if dist < threshold:
                self.highlighted_lines.pop(i)
                key = (x1, y1, x2, y2)
                if key in self.measured_lines: del self.measured_lines[key]
                self.update_overlay()
                return True
        return False

    def highlight_connected_component(self, x, y):
        if self.check_click_on_highlighted_component(x, y): return
        
        if self.current_threshold_binary is not None and self.current_threshold_binary[y, x] > 0:
            num, labels = cv2.connectedComponents(self.current_threshold_binary)
            lbl = labels[y, x]
            mask = (labels == lbl).astype(np.uint8) * 255
            
            for i, (m, _, _) in enumerate(self.highlighted_components):
                if np.array_equal(m, mask):
                    self.highlighted_components.pop(i)
                    self.update_overlay()
                    return

            self.highlighted_components.append((mask, self.threshold_value, self.color_index))
            cid = hash(mask.tobytes())
            _, l_um = self.measure_component_length(mask)
            self.measured_components[cid] = l_um
            print(f"Component measured: {l_um:.2f} um")
            self.update_overlay()
        else:
            print("No component at this location.")

    def get_line_pixels(self, x1, y1, x2, y2):
        pixels = []
        dx = abs(x2 - x1); dy = abs(y2 - y1)
        sx = 1 if x1 < x2 else -1; sy = 1 if y1 < y2 else -1
        err = dx - dy
        x, y = x1, y1
        while True:
            pixels.append((x, y))
            if x == x2 and y == y2: break
            e2 = 2 * err
            if e2 > -dy: err -= dy; x += sx
            if e2 < dx: err += dx; y += sy
        return pixels

    def highlight_line(self, x1, y1, x2, y2):
        self.highlighted_lines.append((x1, y1, x2, y2, self.color_index))
        dist = np.sqrt((x2-x1)**2 + (y2-y1)**2)
        l_um = dist / THRESHOLD_PIXELS_PER_MICRON
        self.measured_lines[(x1, y1, x2, y2)] = l_um
        print(f"Line measured: {l_um:.2f} um")
        self.update_overlay()
    
    def highlight_pixel(self, x, y):
        for i, (px, py, _) in enumerate(self.highlighted_pixels):
            if (px, py) == (x, y):
                self.highlighted_pixels.pop(i)
                if (x,y) in self.measured_pixels: del self.measured_pixels[(x,y)]
                self.update_overlay()
                return
        
        self.highlighted_pixels.append((x, y, self.color_index))
        self.measured_pixels[(x,y)] = 1.0 / THRESHOLD_PIXELS_PER_MICRON
        self.update_overlay()
    
    def save_threshold_image(self):
        out = filedialog.asksaveasfilename(defaultextension=".png")
        if out: cv2.imwrite(out, self.current_threshold_binary)

    def setup_windows(self):
        cv2.namedWindow("Sperm Overlay", cv2.WINDOW_NORMAL)
        cv2.namedWindow("Instructions", cv2.WINDOW_NORMAL)
        cv2.resizeWindow("Sperm Overlay", 1200, 1000)
        cv2.moveWindow("Sperm Overlay", 10, 50)
        cv2.moveWindow("Instructions", 850, 50)
        cv2.resizeWindow("Instructions", 1000, 1000)


    def run_gui(self):
        self.setup_windows()
        cv2.setMouseCallback("Sperm Overlay", self.click_event)

        while self.running:
            if self.threshold_mode:
                display = self.generate_threshold_preview()
                if self.selecting and self.selection_start:
                      d_rect = display.copy()
                      cv2.rectangle(d_rect, self.selection_start, self.selection_end, (255,255,255), 2)
                      cv2.imshow("Sperm Overlay", d_rect)
                else:
                    cv2.imshow("Sperm Overlay", display)
            else:
                # In Auto Mode, show the auto-overlay (or plain overlay if failed)
                if self.auto_overlay is not None:
                    to_show = self.auto_overlay.copy()
                else:
                    to_show = cv2.cvtColor(self.img_gray.copy(), cv2.COLOR_GRAY2BGR)
                
                # Show manual highlights on top of Auto Mode too if enabled
                if self.show_highlights:
                    self.apply_highlights_to_image(to_show)

                cv2.imshow("Sperm Overlay", to_show)
            
            cv2.imshow("Instructions", self.create_instructions_image())
            
            key = cv2.waitKey(1) & 0xFF

            if key == 27: self.running = False
            elif key == ord('i'): 
                self.toggle_inversion()
                print(f"Inversion: {self.inverted}")
            elif key == ord('l'): self.load_new_image()
            elif key == ord('h'): 
                self.show_highlights = not self.show_highlights
                self.update_overlay()
            elif key == ord('m'):
                if not self.threshold_mode:
                    self.run_auto_measure()
            elif key == ord('t'):
                self.threshold_mode = not self.threshold_mode
                print(f"Switched to {'Manual Threshold' if self.threshold_mode else 'Auto Measure'} Mode")
                if not self.threshold_mode: self.update_overlay()
            elif key == ord('p'): self.print_color_sums()
            elif key == ord('e'): self.export_color_sums()
            
            # Auto Mode Controls
            if not self.threshold_mode:
                if key == ord('[') or key == ord('{'):
                    self.auto_threshold_adjust = max(0.1, self.auto_threshold_adjust - 0.1)
                    self.run_auto_measure()
                elif key == ord(']') or key == ord('}'):
                    self.auto_threshold_adjust = min(5.0, self.auto_threshold_adjust + 0.1)
                    self.run_auto_measure()
                
                # New Auto Tools
                elif key == ord('n'):
                    self.auto_add_line_mode = not self.auto_add_line_mode
                    self.auto_add_pixel_mode = False
                    self.auto_remove_mode = False
                    self.auto_line_start = None
                    self.draw_auto_overlay()
                elif key == ord('b'):
                    self.auto_add_pixel_mode = not self.auto_add_pixel_mode
                    self.auto_add_line_mode = False
                    self.auto_remove_mode = False
                    self.auto_line_start = None
                    self.draw_auto_overlay()
                elif key == ord('r'):
                    self.auto_remove_mode = not self.auto_remove_mode
                    self.auto_add_line_mode = False
                    self.auto_add_pixel_mode = False
                    self.auto_line_start = None
                    self.draw_auto_overlay()
            
            # Manual Threshold Controls (Only active in Threshold Mode)
            if self.threshold_mode:
                if key == ord('+') or key == ord('=') or key == 82: 
                    self.threshold_value = min(255, self.threshold_value + 5)
                elif key == ord('-') or key == ord('_') or key == 84: 
                    self.threshold_value = max(0, self.threshold_value - 5)
                elif key == ord('c'):
                    self.highlighted_components = []
                    self.measured_components = {}
                    self.highlighted_pixels = []
                    self.measured_pixels = {}
                    self.highlighted_lines = []
                    self.measured_lines = {}
                    print("Cleared all.")
                elif key == ord('u'):
                    self.pixel_measuring_mode = not self.pixel_measuring_mode
                    print(f"Pixel Mode: {self.pixel_measuring_mode}")
                elif key == ord('k'):
                    self.line_mode = not self.line_mode
                    print(f"Line Tool Mode: {self.line_mode}")
                elif key == ord('s'): self.save_threshold_image()
                elif ord('0') <= key <= ord('9'):
                    digit = int(chr(key))
                    self.color_index = 9 if digit == 0 else digit - 1
                    print(f"Color: {self.color_index+1}")

        cv2.destroyAllWindows()
        self.export_results()

    def smooth_line(self, binary):
        blurred = cv2.GaussianBlur(binary, (5, 5), 0)
        _, binary_smooth = cv2.threshold(blurred, 0, 255, cv2.THRESH_OTSU)
        return binary_smooth

    def skeletonize_cell(self, binary):
        skel = skeletonize(binary > 0).astype(np.uint8)
        labels, counts = np.unique(cv2.connectedComponents(skel)[1], return_counts=True)
        if len(counts) > 1:
            largest = labels[1:][np.argmax(counts[1:])] 
            skel = (cv2.connectedComponents(skel)[1] == largest).astype(np.uint8)
        return skel

    def trim_endpoints(self, skeleton, points, trim_radius=None):
        if trim_radius is None: trim_radius = THRESHOLD_TRIM_RADIUS
        ys, xs = np.nonzero(skeleton)
        coords = np.stack([xs, ys], axis=1)
        for point in [points[0], points[-1]]:
            dists = cdist(coords, np.array([point]))
            to_trim = dists < trim_radius
            coords = coords[~to_trim.ravel()]
        new_skeleton = np.zeros_like(skeleton)
        for x, y in coords: new_skeleton[int(y), int(x)] = 1
        return new_skeleton

    def measure_length(self, skeleton, pixels_per_micron=None):
        if pixels_per_micron is None: pixels_per_micron = THRESHOLD_PIXELS_PER_MICRON
        ys, xs = np.nonzero(skeleton)
        coords = np.stack([xs, ys], axis=1)
        if len(coords) < 5: return 0, 0
        try:
            from scipy.interpolate import splprep, splev
            tck, u = splprep(coords.T, s=5)
            x_s, y_s = splev(np.linspace(0,1,400), tck)
            p = np.stack([x_s, y_s], axis=1)
            d = np.sqrt(np.sum(np.diff(p, axis=0)**2, axis=1)).sum()
            return d, d / pixels_per_micron
        except:
             return 0, 0

    def export_results(self, output_csv="results.csv"):
        rows = []
        for cid, length_um in self.measured_components.items():
             rows.append({'Type': 'Component', 'Length_um': length_um})
        for lid, length_um in self.measured_lines.items():
             rows.append({'Type': 'Line', 'Length_um': length_um})
        
        if rows:
            df = pd.DataFrame(rows)
            df.to_csv(output_csv, index=False)
            print(f"Results saved to {output_csv}")

if __name__ == "__main__":
    Tk().withdraw()
    print("Select an image file...")
    img_path = filedialog.askopenfilename(
        title="Choose sperm image",
        filetypes=[("Images", "*.jpg *.jpeg *.png *.bmp *.tif *.tiff"), ("All files", "*")]
    )

    if not img_path:
        print("No file selected. Exiting.")
        exit()

    print(f"Selected: {img_path}")
    sizer = SpermSizer(img_path)
    sizer.run_gui()